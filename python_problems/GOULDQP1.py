from s2xlib import *
class  GOULDQP1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GOULDQP1
#    *********
# 
#    Source: problem 118 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981, as modified by N.I.M. Gould in "An algorithm 
#    for large-scale quadratic programming", IMA J. Num. Anal (1991),
#    11, 299-324, problem class 1.
# 
#    SIF input: B Baudson, Jan 1990 modified by Nick Gould, Jan, 2011
# 
#    classification = "QLR2-AN-32-17"
# 
#    Other useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GOULDQP1'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'GOULDQP1'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['4'] = 4
        v_['5'] = 5
        v_['8'] = 8
        v_['12'] = 12
        v_['15'] = 15
        v_['17'] = 17
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['15'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for K in range(int(v_['1']),int(v_['4'])+1):
            [iv,ix_,_] = s2x_ii('AS'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'AS'+str(K))
        for K in range(int(v_['1']),int(v_['4'])+1):
            [iv,ix_,_] = s2x_ii('CS'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'CS'+str(K))
        for K in range(int(v_['1']),int(v_['4'])+1):
            [iv,ix_,_] = s2x_ii('BS'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'BS'+str(K))
        for K in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2x_ii('DS'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'DS'+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for K in range(int(v_['0']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['3K+1']))]
            pbm.A[ig,iv] = float(2.3)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K+2']))]
            pbm.A[ig,iv] = float(1.7)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K+3']))]
            pbm.A[ig,iv] = float(2.2)+pbm.A[ig,iv]
        for K in range(int(v_['1']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            v_['3K-2'] = -2+v_['3K']
            v_['3K-1'] = -1+v_['3K']
            [ig,ig_,_] = s2x_ii('A'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(K))
            iv = ix_['X'+str(int(v_['3K+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K-2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['AS'+str(K)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('B'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'B'+str(K))
            iv = ix_['X'+str(int(v_['3K+3']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['BS'+str(K)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(K))
            iv = ix_['X'+str(int(v_['3K+2']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['CS'+str(K)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('D1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DS1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('D2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D2')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DS2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('D3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D3')
        iv = ix_['X7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DS3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('D4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D4')
        iv = ix_['X10']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DS4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('D5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D5')
        iv = ix_['X13']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['DS5']
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
        for K in range(int(v_['1']),int(v_['4'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['A'+str(K)],float(-7.0))
            pbm.gconst = arrset(pbm.gconst,ig_['B'+str(K)],float(-7.0))
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(K)],float(-7.0))
        pbm.gconst = arrset(pbm.gconst,ig_['D1'],float(60.0))
        pbm.gconst = arrset(pbm.gconst,ig_['D2'],float(50.0))
        pbm.gconst = arrset(pbm.gconst,ig_['D3'],float(70.0))
        pbm.gconst = arrset(pbm.gconst,ig_['D4'],float(85.0))
        pbm.gconst = arrset(pbm.gconst,ig_['D5'],float(100.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 8.0
        pb.xupper[ix_['X1']] = 21.0
        pb.xlower[ix_['X2']] = 43.0
        pb.xupper[ix_['X2']] = 57.0
        pb.xlower[ix_['X3']] = 3.0
        pb.xupper[ix_['X3']] = 16.0
        for K in range(int(v_['1']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            pb.xupper[ix_['X'+str(int(v_['3K+1']))]] = 90.0
            pb.xupper[ix_['X'+str(int(v_['3K+2']))]] = 120.0
            pb.xupper[ix_['X'+str(int(v_['3K+3']))]] = 60.0
        for K in range(int(v_['1']),int(v_['4'])+1):
            pb.xlower[ix_['AS'+str(K)]] = 0.0
            pb.xupper[ix_['AS'+str(K)]] = 13.0
            pb.xlower[ix_['BS'+str(K)]] = 0.0
            pb.xupper[ix_['BS'+str(K)]] = 13.0
            pb.xlower[ix_['CS'+str(K)]] = 0.0
            pb.xupper[ix_['CS'+str(K)]] = 14.0
        pb.xlower[ix_['DS1']] = 0.0
        pb.xlower[ix_['DS2']] = 0.0
        pb.xlower[ix_['DS3']] = 0.0
        pb.xlower[ix_['DS4']] = 0.0
        pb.xlower[ix_['DS5']] = 0.0
        pb.xupper[ix_['DS1']] = 60.0
        pb.xupper[ix_['DS2']] = 50.0
        pb.xupper[ix_['DS3']] = 70.0
        pb.xupper[ix_['DS4']] = 85.0
        pb.xupper[ix_['DS5']] = 100.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(20.0))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(55.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(55.0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(15.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(15.0)))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(60.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(60.0)))
        if('X8' in ix_):
            pb.x0[ix_['X8']] = float(60.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8']),float(60.0)))
        if('X11' in ix_):
            pb.x0[ix_['X11']] = float(60.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X11']),float(60.0)))
        if('X14' in ix_):
            pb.x0[ix_['X14']] = float(60.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X14']),float(60.0)))
        if('AS1' in ix_):
            pb.x0[ix_['AS1']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['AS1']),float(7.0)))
        if('BS1' in ix_):
            pb.x0[ix_['BS1']] = float(12.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['BS1']),float(12.0)))
        if('CS1' in ix_):
            pb.x0[ix_['CS1']] = float(12.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['CS1']),float(12.0)))
        if('AS2' in ix_):
            pb.x0[ix_['AS2']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['AS2']),float(7.0)))
        if('BS2' in ix_):
            pb.x0[ix_['BS2']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['BS2']),float(7.0)))
        if('CS2' in ix_):
            pb.x0[ix_['CS2']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['CS2']),float(7.0)))
        if('AS3' in ix_):
            pb.x0[ix_['AS3']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['AS3']),float(7.0)))
        if('BS3' in ix_):
            pb.x0[ix_['BS3']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['BS3']),float(7.0)))
        if('CS3' in ix_):
            pb.x0[ix_['CS3']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['CS3']),float(7.0)))
        if('AS4' in ix_):
            pb.x0[ix_['AS4']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['AS4']),float(7.0)))
        if('BS4' in ix_):
            pb.x0[ix_['BS4']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['BS4']),float(7.0)))
        if('CS4' in ix_):
            pb.x0[ix_['CS4']] = float(7.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['CS4']),float(7.0)))
        if('DS1' in ix_):
            pb.x0[ix_['DS1']] = float(30.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['DS1']),float(30.0)))
        if('DS2' in ix_):
            pb.x0[ix_['DS2']] = float(50.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['DS2']),float(50.0)))
        if('DS3' in ix_):
            pb.x0[ix_['DS3']] = float(30.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['DS3']),float(30.0)))
        if('DS4' in ix_):
            pb.x0[ix_['DS4']] = float(15.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['DS4']),float(15.0)))
        if('DS5' in ix_):
            pb.x0[ix_['DS5']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['DS5']),float(0.0)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['15'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,20.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
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
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.00015))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(10.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(25.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.5))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.00015))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.0001))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.00015))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "QLR2-AN-32-17"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

