from s2xlib import *
class  HS84(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    *********
#    Problem : HS84
#    *********
# 
#    Source: problem 84 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "QQR2-AN-5-3"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS84'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS84'
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
        v_['15'] = 15
        v_['16'] = 16
        v_['17'] = 17
        v_['18'] = 18
        v_['19'] = 19
        v_['20'] = 20
        v_['21'] = 21
        v_['A'+str(int(v_['1']))] = -24345.0
        v_['A'+str(int(v_['2']))] = -8720288.849
        v_['MA'+str(int(v_['2']))] = -1.0*v_['A'+str(int(v_['2']))]
        v_['A'+str(int(v_['3']))] = 150512.5253
        v_['MA'+str(int(v_['3']))] = -1.0*v_['A'+str(int(v_['3']))]
        v_['A'+str(int(v_['4']))] = -156.6950325
        v_['MA'+str(int(v_['4']))] = -1.0*v_['A'+str(int(v_['4']))]
        v_['A'+str(int(v_['5']))] = 476470.3222
        v_['MA'+str(int(v_['5']))] = -1.0*v_['A'+str(int(v_['5']))]
        v_['A'+str(int(v_['6']))] = 729482.8271
        v_['MA'+str(int(v_['6']))] = -1.0*v_['A'+str(int(v_['6']))]
        v_['A'+str(int(v_['7']))] = -145421.402
        v_['A'+str(int(v_['8']))] = 2931.1506
        v_['A'+str(int(v_['9']))] = -40.427932
        v_['A'+str(int(v_['10']))] = 5106.192
        v_['A'+str(int(v_['11']))] = 15711.36
        v_['A'+str(int(v_['12']))] = -155011.1084
        v_['A'+str(int(v_['13']))] = 4360.53352
        v_['A'+str(int(v_['14']))] = 12.9492344
        v_['A'+str(int(v_['15']))] = 10236.884
        v_['A'+str(int(v_['16']))] = 13176.786
        v_['A'+str(int(v_['17']))] = -326669.5104
        v_['A'+str(int(v_['18']))] = 7390.68412
        v_['A'+str(int(v_['19']))] = -27.8986976
        v_['A'+str(int(v_['20']))] = 16643.076
        v_['A'+str(int(v_['21']))] = 30988.146
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
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
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['MA'+str(int(v_['2']))])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['A'+str(int(v_['7']))])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON2')
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['A'+str(int(v_['12']))])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON3')
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['A'+str(int(v_['17']))])+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(v_['A'+str(int(v_['1']))]))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        grange = arrset(grange,ig_['CON1'],float(294000.0))
        grange = arrset(grange,ig_['CON2'],float(294000.0))
        grange = arrset(grange,ig_['CON3'],float(277200.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X'+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['2']))]] = 1.2
        pb.xlower[ix_['X'+str(int(v_['3']))]] = 20.0
        pb.xlower[ix_['X'+str(int(v_['4']))]] = 9.0
        pb.xlower[ix_['X'+str(int(v_['5']))]] = 6.5
        pb.xupper[ix_['X'+str(int(v_['1']))]] = 1000.0
        pb.xupper[ix_['X'+str(int(v_['2']))]] = 2.4
        pb.xupper[ix_['X'+str(int(v_['3']))]] = 60.0
        pb.xupper[ix_['X'+str(int(v_['4']))]] = 9.3
        pb.xupper[ix_['X'+str(int(v_['5']))]] = 7.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['1']))]] = float(2.52)
        pb.x0[ix_['X'+str(int(v_['2']))]] = float(2.0)
        pb.x0[ix_['X'+str(int(v_['3']))]] = float(37.5)
        pb.x0[ix_['X'+str(int(v_['4']))]] = float(9.25)
        pb.x0[ix_['X'+str(int(v_['5']))]] = float(6.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['4'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            v_['IP1'] = 1+I
            vname = 'X'+str(int(v_['1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['IP1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MA'+str(int(v_['3']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MA'+str(int(v_['4']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MA'+str(int(v_['5']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MA'+str(int(v_['6']))]))
        ig = ig_['CON1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['8']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['9']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['10']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['11']))]))
        ig = ig_['CON2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['13']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['14']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['15']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['16']))]))
        ig = ig_['CON3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['18']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['19']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['20']))]))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['21']))]))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        pb.cupper[np.arange(pb.nge)] = grange[gegrps]
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QQR2-AN-5-3"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

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
