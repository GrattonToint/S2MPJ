from s2xlib import *
class  HS67(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS67
#    *********
# 
#    Source: problem 67 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 8 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn & Nick Gould, April 1991.
# 
#    classification = "OOI2-AN-3-14"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS67'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS67'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['A'+str(v_['1'])] = 0.0
        v_['A'+str(v_['2'])] = 0.0
        v_['A'+str(v_['3'])] = 85.0
        v_['A'+str(v_['4'])] = 90.0
        v_['A'+str(v_['5'])] = 3.0
        v_['A'+str(v_['6'])] = 0.01
        v_['A'+str(v_['7'])] = 145.0
        v_['A'+str(v_['8'])] = 5000.0
        v_['A'+str(v_['9'])] = 2000.0
        v_['A'+str(v_['10'])] = 93.0
        v_['A'+str(v_['11'])] = 95.0
        v_['A'+str(v_['12'])] = 12.0
        v_['A'+str(v_['13'])] = 4.0
        v_['A'+str(v_['14'])] = 162.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['3'])+1):
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
        iv = ix_['X'+str(v_['1'])]
        pbm.A[ig,iv] = 5.04+pbm.A[ig,iv]
        iv = ix_['X'+str(v_['2'])]
        pbm.A[ig,iv] = 0.035+pbm.A[ig,iv]
        iv = ix_['X'+str(v_['3'])]
        pbm.A[ig,iv] = 10.0+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['7'])+1):
            [ig,ig_,_] = s2x_ii('AG'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'AG'+str(I))
            [ig,ig_,_] = s2x_ii('AL'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'AL'+str(I))
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
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+7'] = 7+I
            pbm.gconst = arrset(pbm.gconst,ig_['AG'+str(I)],v_['A'+str(I)])
            pbm.gconst = arrset(pbm.gconst,ig_['AL'+str(I)],v_['A'+str(v_['I+7'])])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.00001)
        pb.xupper = np.full((pb.n,1),+float('inf'))
        pb.xupper[ix_['X'+str(v_['1'])]] = 2000.0
        pb.xupper[ix_['X'+str(v_['2'])]] = 16000.0
        pb.xupper[ix_['X'+str(v_['3'])]] = 120.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(v_['1'])]] = 1745.0
        pb.x0[ix_['X'+str(v_['2'])]] = 12000.0
        pb.x0[ix_['X'+str(v_['3'])]] = 110.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'Y2Y5', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y2', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y3', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y4', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y5', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y6', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y7', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2x_ii( 'Y8', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E25'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y2Y5')
        ielftype = arrset(ielftype, ie, iet_["Y2Y5"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y2')
        ielftype = arrset(ielftype, ie, iet_["Y2"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y3')
        ielftype = arrset(ielftype, ie, iet_["Y3"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y4')
        ielftype = arrset(ielftype, ie, iet_["Y4"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y5')
        ielftype = arrset(ielftype, ie, iet_["Y5"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y6')
        ielftype = arrset(ielftype, ie, iet_["Y6"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y7')
        ielftype = arrset(ielftype, ie, iet_["Y7"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'Y8')
        ielftype = arrset(ielftype, ie, iet_["Y8"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.00001,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,-0.063)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,3.36)
        ig = ig_['AG'+str(v_['1'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['1'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['2'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['2'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['3'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['3'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['4'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['4'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['5'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['5'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['6'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['6'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AG'+str(v_['7'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AL'+str(v_['7'])]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOI2-AN-3-14"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
