from s2xlib import *
class  HS85(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS85
#    *********
# 
#    The problem is to optimize the net profit of an hypothetical
#    wood-pulp plant. The constraints include the usual material
#    and energy balances as well as several empirical equations.
# 
#    Source: problem 85 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#      SAVEs removed December 3rd 2014
# 
#    classification = "OOI2-MN-5-21"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS85'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS85'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 5
        v_['A2'] = 17.505
        v_['A3'] = 11.275
        v_['A4'] = 214.228
        v_['A5'] = 7.458
        v_['A6'] = 0.961
        v_['A7'] = 1.612
        v_['A8'] = 0.146
        v_['A9'] = 107.99
        v_['A10'] = 922.693
        v_['A11'] = 926.832
        v_['A12'] = 18.766
        v_['A13'] = 1072.163
        v_['A14'] = 8961.448
        v_['A15'] = 0.063
        v_['A16'] = 71084.33
        v_['A17'] = 2802713.0
        v_['B2'] = 1053.6667
        v_['B3'] = 35.03
        v_['B4'] = 665.585
        v_['B5'] = 584.463
        v_['B6'] = 265.916
        v_['B7'] = 7.046
        v_['B8'] = 0.222
        v_['B9'] = 273.366
        v_['B10'] = 1286.105
        v_['B11'] = 1444.046
        v_['B12'] = 537.141
        v_['B13'] = 3247.039
        v_['B14'] = 26844.086
        v_['B15'] = 0.386
        v_['B16'] = 140000.0
        v_['B17'] = 12146108.0
        v_['1'] = 1
        v_['2'] = 2
        v_['17'] = 17
        v_['19'] = 19
        v_['20'] = 20
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
        [ig,ig_,_] = s2x_ii('CON0',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON0')
        iv = ix_['X2']
        pbm.A[ig,iv] = 1.5+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = -1.0+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['20'])+1):
            [ig,ig_,_] = s2x_ii('CON'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CON'+str(I))
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],0.1365)
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],213.1)
        for I in range(int(v_['2']),int(v_['17'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CON'+str(I)],v_['A'+str(I)])
        pbm.gconst = arrset(pbm.gconst,ig_['CON19'],-21.0)
        pbm.gconst = arrset(pbm.gconst,ig_['CON20'],110.6)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        for I in range(int(v_['2']),int(v_['17'])+1):
            v_['DIF'] = v_['B'+str(I)]-v_['A'+str(I)]
            grange = arrset(grange,ig_['CON'+str(I)],v_['DIF'])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 7.044148e+2
        pb.xupper[ix_['X1']] = 9.063855e+2
        pb.xlower[ix_['X2']] = 6.86e+1
        pb.xupper[ix_['X2']] = 2.8888e+2
        pb.xupper[ix_['X3']] = 1.3475e+2
        pb.xlower[ix_['X4']] = 1.930e+2
        pb.xupper[ix_['X4']] = 2.870966e+2
        pb.xlower[ix_['X5']] = 2.50e+1
        pb.xupper[ix_['X5']] = 8.41988e+1
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = 9.0e+2
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),9.0e+2)
        if('X2' in ix_):
            pb.x0[ix_['X2']] = 8.0e+1
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),8.0e+1)
        if('X3' in ix_):
            pb.x0[ix_['X3']] = 1.15e+2
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),1.15e+2)
        if('X4' in ix_):
            pb.x0[ix_['X4']] = 2.67e+2
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),2.67e+2)
        if('X5' in ix_):
            pb.x0[ix_['X5']] = 2.7e+1
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),2.7e+1)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'Y', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftp = []
        elftp = loaset(elftp,it,0,'PI')
        [it,iet_,_] = s2x_ii( 'C', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftp = loaset(elftp,it,0,'PI')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['19'])+1):
            v_['PI'] = I
            v_['PI'] = 0.01+v_['PI']
            ename = 'C'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'C')
            ielftype = arrset(ielftype, ie, iet_["C"])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='PI')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['PI'])
        for I in range(int(v_['1']),int(v_['20'])+1):
            v_['PI'] = I
            v_['PI'] = 0.01+v_['PI']
            ename = 'Y'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'Y')
            ielftype = arrset(ielftype, ie, iet_["Y"])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='PI')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['PI'])
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y17'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,-5.843e-7)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y14'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.17e-4)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,2.358e-5)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y16'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.502e-6)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,0.0321)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,0.00423)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C18'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.0e-4)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C19'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,37.48)
        for I in range(int(v_['1']),int(v_['20'])+1):
            ig = ig_['CON'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "OOI2-MN-5-21"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
