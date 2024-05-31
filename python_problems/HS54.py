from s2xlib import *
class  HS54(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS54
#    *********
# 
#    Source: problem 54, incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Betts problem 11.7, JOTA 21, 1977, pp.137-174.
#    SIF input: A.R. Conn, April 1990 and Nick Gould, October 1990
# 
#    classification = "OLR2-AN-6-1"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS54'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS54'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 6
        v_['1'] = 1
        v_['6'] = 6
        v_['RHO'] = 2.0e-1
        v_['RHOSQR'] = v_['RHO']*v_['RHO']
        v_['1-RHOSQR'] = 1.0-v_['RHOSQR']
        v_['FACTOR'] = 1.0/v_['1-RHOSQR']
        v_['MU1'] = 1.0e+4
        v_['MU2'] = 1.0e+0
        v_['MU3'] = 2.0e+6
        v_['MU4'] = 1.0e+1
        v_['MU5'] = 1.0e-3
        v_['MU6'] = 1.0e+8
        v_['SIGMA1'] = 8.0e+3
        v_['SIGMA2'] = 1.0e+0
        v_['SIGMA3'] = 7.0e+6
        v_['SIGMA4'] = 5.0e+1
        v_['SIGMA5'] = 5.0e-2
        v_['SIGMA6'] = 5.0e+8
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
        [ig,ig_,_] = s2x_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(4.0e+3)+pbm.A[ig,iv]
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
        v_['0.2SI1'] = 2.0e-1*v_['SIGMA1']
        v_['2000SI2'] = 2.0e+3*v_['SIGMA2']
        v_['4000MU2'] = 4.0e+3*v_['MU2']
        v_['RHS'] = v_['MU1']+v_['4000MU2']
        v_['RHS'] = v_['RHS']+v_['0.2SI1']
        v_['RHS'] = v_['RHS']+v_['2000SI2']
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(v_['RHS']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['X1']] = 2.0e+4
        pb.xlower[ix_['X2']] = - 1.0e+1
        pb.xupper[ix_['X2']] = 1.0e+1
        pb.xupper[ix_['X3']] = 1.0e+7
        pb.xupper[ix_['X4']] = 2.0e+1
        pb.xlower[ix_['X5']] = - 1.0e+0
        pb.xupper[ix_['X5']] = 1.0e+0
        pb.xupper[ix_['X6']] = 2.0e+8
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(6.0e+3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(6.0e+3)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(1.5e+0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(1.5e+0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(4.0e+6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(4.0e+6)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(2.0e+0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4']),float(2.0e+0)))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(3.0e-3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(3.0e-3)))
        if('X6' in ix_):
            pb.x0[ix_['X6']] = float(5.0e+7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X6']),float(5.0e+7)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftp = []
        elftp = loaset(elftp,it,0,'MU')
        elftp = loaset(elftp,it,1,'SIGMA')
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'MU1')
        elftp = loaset(elftp,it,1,'MU2')
        elftp = loaset(elftp,it,2,'SIGMA1')
        elftp = loaset(elftp,it,3,'SIGMA2')
        elftp = loaset(elftp,it,4,'RHO')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['6'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
                ielftype = arrset( ielftype,ie,iet_['eSQR'])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='MU')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['MU'+str(I)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='SIGMA')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['SIGMA'+str(I)]))
        ename = 'F1'
        [ie,ie_,newelt] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='RHO')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RHO']))
        posep = find(elftp[ielftype[ie]],lambda x:x=='MU1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['MU1']))
        posep = find(elftp[ielftype[ie]],lambda x:x=='MU2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['MU2']))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SIGMA1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['SIGMA1']))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SIGMA2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['SIGMA2']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gNORMAL',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        pbm.grftype = arrset(pbm.grftype,ig,'gNORMAL')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['FACTOR']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['FACTOR']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0e+0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0e+0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0e+0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0e+0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['FACTOR']))
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
        pb.pbclass = "OLR2-AN-6-1"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V1MP = (EV_[0]-pbm.elpar[iel_][0])/pbm.elpar[iel_][1]
        f_   = V1MP**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*V1MP/pbm.elpar[iel_][1]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0/pbm.elpar[iel_][1]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TERM1 = (EV_[0]-pbm.elpar[iel_][0])/pbm.elpar[iel_][2]
        TERM2 = (EV_[1]-pbm.elpar[iel_][1])/pbm.elpar[iel_][3]
        RHO2 = pbm.elpar[iel_][4]+pbm.elpar[iel_][4]
        f_   = RHO2*TERM1*TERM2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RHO2*TERM2/pbm.elpar[iel_][2]
            g_[1] = RHO2*TERM1/pbm.elpar[iel_][3]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = RHO2/(pbm.elpar[iel_][2]*pbm.elpar[iel_][3])
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gNORMAL(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        EXPHV = np.exp(-0.5*GVAR_)
        f_= -EXPHV
        if nargout>1:
            g_ = 5.0e-1*EXPHV
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -2.5e-1*EXPHV
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

