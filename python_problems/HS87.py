from s2xlib import *
class  HS87(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS87
#    *********
# 
#    Optimization of an electrical network (EDF) by P. Huard.
# 
#    Source: problem 87 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Note: There are two variants described in the papers
# 
#       D.H. Himmelblau "Applied nonlinear programming",
#       McGraw-Hill, New-York, 1972, problem 15,
# 
#    and
# 
#       A.R. Colville, "A comparative study on nonlinear programming",
#       IBM Scientific Center Report 320-2949, New York, 1968, problem 6.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "OOI2-MN-6-4"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS87'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS87'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['A'] = 131.078
        v_['B'] = 1.48577
        v_['C'] = 0.90798
        v_['D'] = np.cos(1.47588)
        v_['E'] = np.sin(1.47588)
        v_['F'] = 1.48577
        v_['1'] = 1
        v_['-B'] = -1.0e+0*v_['B']
        v_['-F'] = -1.0e+0*v_['F']
        v_['C/A'] = v_['C']/v_['A']
        v_['1/A'] = 1.0e+0/v_['A']
        v_['-1/A'] = -1.0e+0*v_['1/A']
        v_['CD/A'] = v_['C/A']*v_['D']
        v_['CE/A'] = v_['C/A']*v_['E']
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
        [ig,ig_,_] = s2x_ii('C1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C2')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C3')
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('C4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'C4')
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(-3.0e+2))
        pbm.gconst = arrset(pbm.gconst,ig_['C4'],float(-2.0e+2))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X3']] = 340.0
        pb.xlower[ix_['X4']] = 340.0
        pb.xlower[ix_['X5']] = -1000.0
        pb.xupper[ix_['X1']] = 400.0
        pb.xupper[ix_['X2']] = 1000.0
        pb.xupper[ix_['X3']] = 420.0
        pb.xupper[ix_['X4']] = 420.0
        pb.xupper[ix_['X5']] = 10000.0
        pb.xupper[ix_['X6']] = 0.5236
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(107.8119)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(107.8119)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(196.3186)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),float(196.3186)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(373.8307)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(373.8307)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(420.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),float(420.0))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(21.30713)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(21.30713)))
        if('X6' in ix_):
            pb.x0[ix_['X6']] = float(0.153292)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X6']),float(0.153292)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eF1', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2x_ii( 'eF2', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2x_ii( 'eCOS', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2x_ii( 'eSIN', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2x_ii( 'eSQUARE', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'OF1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eF1')
        ielftype = arrset(ielftype, ie, iet_["eF1"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'OF2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eF2')
        ielftype = arrset(ielftype, ie, iet_["eF2"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C1E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS')
        ielftype = arrset(ielftype, ie, iet_["eCOS"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-F']))
        ename = 'C1E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eSQUARE"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C2E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOS')
        ielftype = arrset(ielftype, ie, iet_["eCOS"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
        ename = 'C2E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eSQUARE"])
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C3E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN')
        ielftype = arrset(ielftype, ie, iet_["eSIN"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
        ename = 'C3E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eSQUARE"])
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'C4E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSIN')
        ielftype = arrset(ielftype, ie, iet_["eSIN"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='P')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-B']))
        ename = 'C4E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQUARE')
        ielftype = arrset(ielftype, ie, iet_["eSQUARE"])
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OF1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OF2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['C1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C1E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/A']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C1E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CD/A']))
        ig = ig_['C2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C2E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/A']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C2E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CD/A']))
        ig = ig_['C3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C3E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/A']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C3E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CE/A']))
        ig = ig_['C4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C4E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/A']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C4E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CE/A']))
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
        pb.pbclass = "OOI2-MN-6-4"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eF1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        I1 = EV_[0]<300.0
        I2 = EV_[0]>=300.0
        if I1!=0:
            F = 30.0*EV_[0]
        if I2!=0:
            F = 31.0*EV_[0]
        if I1!=0:
            G = 30.0
        if I2!=0:
            G = 31.0
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eF2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        I1 = EV_[0]<100.0
        I2 = EV_[0]>=100.0 and EV_[0]<200.0
        I3 = EV_[0]>=200.0
        if I1!=0:
            F = 28.0*EV_[0]
        if I2!=0:
            F = 29.0*EV_[0]
        if I3!=0:
            F = 30.0*EV_[0]
        if I1!=0:
            G = 28.0
        if I2!=0:
            G = 29.0
        if I3!=0:
            G = 30.0
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SN = np.sin(EV_[2]+pbm.elpar[iel_][0])
        CS = np.cos(EV_[2]+pbm.elpar[iel_][0])
        f_   = EV_[0]*EV_[1]*CS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*CS
            g_[1] = EV_[0]*CS
            g_[2] = -EV_[0]*EV_[1]*SN
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = CS
                H_[1,0] = H_[0,1]
                H_[0,2] = -EV_[1]*SN
                H_[2,0] = H_[0,2]
                H_[1,2] = -EV_[0]*SN
                H_[2,1] = H_[1,2]
                H_[2,2] = -EV_[0]*EV_[1]*CS
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSIN(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SN = np.sin(EV_[2]+pbm.elpar[iel_][0])
        CS = np.cos(EV_[2]+pbm.elpar[iel_][0])
        f_   = EV_[0]*EV_[1]*SN
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*SN
            g_[1] = EV_[0]*SN
            g_[2] = EV_[0]*EV_[1]*CS
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = SN
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*CS
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*CS
                H_[2,1] = H_[1,2]
                H_[2,2] = -EV_[0]*EV_[1]*SN
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQUARE(pbm,nargout,*args):

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
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

