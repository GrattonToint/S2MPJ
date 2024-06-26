from s2xlib import *
class  ORTHRDM2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHRDM2
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) a planar curve to a set of points
#    in the plane. This set of points is generated by perturbing a
#    first set lying exactly on the predefined curve.
#    The curve is referred to as a cardioid in the original paper,
#    but is in fact a circle.
# 
#    This problem is a modification of ORTHREGD.SIF.
#    The start point was moved to keep any data points from
#    converging to the center of the circle, because such
#    points make the constraint Jacobian singular.
# 
#    Source: adapted from:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, Mar 1991,
#               modified by T. Plantagena, May 1994.
# 
#    classification = "QOR2-AY-V-V"
# 
#    Number of data points
#    (number of variables = 2 NPTS + 3 )
# 
# IE NPTS                100            $-PARAMETER n = 203
# IE NPTS                2000           $-PARAMETER n = 4003     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHRDM2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'ORTHRDM2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NPTS'] = int(4000);  #  SIF file default value
        else:
            v_['NPTS'] = int(args[0])
        v_['TZ3'] = 1.7
        v_['PSEED'] = 237.1531
        v_['PSIZE'] = 0.1
        v_['1'] = 1
        v_['0'] = 0
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['RNPTS'] = float(v_['NPTS'])
        v_['ICR0'] = 1.0/v_['RNPTS']
        v_['INCR'] = v_['ICR0']*v_['2PI']
        v_['Z3SQ'] = v_['TZ3']*v_['TZ3']
        v_['1+TZ3SQ'] = 1.0+v_['Z3SQ']
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['THETA'] = v_['RI-1']*v_['INCR']
            v_['ST'] = np.sin(v_['THETA'])
            v_['CT'] = np.cos(v_['THETA'])
            v_['FACT'] = v_['1+TZ3SQ']+v_['CT']
            v_['R1'] = v_['FACT']*v_['CT']
            v_['R2'] = v_['FACT']*v_['ST']
            v_['XSEED'] = v_['THETA']*v_['PSEED']
            v_['SSEED'] = np.cos(v_['XSEED'])
            v_['PER-1'] = v_['PSIZE']*v_['SSEED']
            v_['PERT'] = 1.0+v_['PER-1']
            v_['XD'+str(I)] = v_['R1']*v_['PERT']
            v_['YD'+str(I)] = v_['R2']*v_['PERT']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('Z1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Z1')
        [iv,ix_,_] = s2x_ii('Z2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Z2')
        [iv,ix_,_] = s2x_ii('Z3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Z3')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2x_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [ig,ig_,_] = s2x_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('OY'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['Y'+str(I)]
            pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['OX'+str(I)],v_['XD'+str(I)])
            pbm.gconst = arrset(pbm.gconst,ig_['OY'+str(I)],v_['YD'+str(I)])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('Z1' in ix_):
            pb.x0[ix_['Z1']] = 1.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Z1']),1.0)
        if('Z2' in ix_):
            pb.x0[ix_['Z2']] = 0.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Z2']),0.0)
        if('Z3' in ix_):
            pb.x0[ix_['Z3']] = 1.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Z3']),1.0)
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            if('X'+str(I) in ix_):
                pb.x0[ix_['X'+str(I)]] = v_['XD'+str(I)]
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)]),v_['XD'+str(I)]))
            if('Y'+str(I) in ix_):
                pb.x0[ix_['Y'+str(I)]] = v_['YD'+str(I)]
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y'+str(I)]),v_['YD'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eTA', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'ZA')
        elftv = loaset(elftv,it,3,'ZB')
        [it,iet_,_] = s2x_ii( 'eTB', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'ZA')
        elftv = loaset(elftv,it,3,'ZB')
        elftv = loaset(elftv,it,4,'ZC')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eTA')
            ielftype = arrset(ielftype, ie, iet_["eTA"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ZA')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ZB')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eTB')
            ielftype = arrset(ielftype, ie, iet_["eTB"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ZA')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ZB')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ZC')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ig = ig_['OX'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['OY'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        pb.pbclass = "QOR2-AY-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eTA(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        T = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        f_   = T*T
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*T*IV_[0]
            g_[1] = 4.0*T*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 4.0*(T+2.0*IV_[0]*IV_[0])
                H_[0,1] = 8.0*IV_[0]*IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 4.0*(T+2.0*IV_[1]*IV_[1])
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,5))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,4] = U_[2,4]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        T = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        ZZSQ = IV_[2]*IV_[2]
        T1 = 1.0+ZZSQ
        T1SQ = T1*T1
        f_   = T*T1SQ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*IV_[0]*T1SQ
            g_[1] = 2.0*IV_[1]*T1SQ
            g_[2] = 4.0*T*T1*IV_[2]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 2.0*T1SQ
                H_[0,2] = 8.0*IV_[0]*T1*IV_[2]
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*T1SQ
                H_[1,2] = 8.0*IV_[1]*T1*IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0*T*(2.0*ZZSQ+T1)
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

