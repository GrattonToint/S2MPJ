from s2xlib import *
class  ORTHREGF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGF
#    *********
# 
#    An orthogonal regression problem
# 
#    The problem is to fit (orthogonally) an torus to a
#    set of points in 3D space. This set of points is generated by
#    perturbing a first set lying exactly on a predefined torus
#    centered at the origin.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "QOR2-AY-V-V"
# 
#    square root of the number of data points
#    (number of variables = 3 * NPTS**2 + 5 )
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHREGF'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'ORTHREGF'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NPTS'] = int(5);  #  SIF file default value
        else:
            v_['NPTS'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE NPTS                7              $-PARAMETER n = 152
# IE NPTS                10             $-PARAMETER n = 305
# IE NPTS                15             $-PARAMETER n = 680
# IE NPTS                20             $-PARAMETER n = 1205
# IE NPTS                40             $-PARAMETER n = 4805
        v_['TP4'] = 1.7
        v_['TP5'] = 0.8
        v_['PSEED'] = 237.1531
        v_['PSIZE'] = 0.2
        v_['1'] = 1
        v_['5'] = 5
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['RNPTS'] = float(v_['NPTS'])
        v_['ICR0'] = 1.0/v_['RNPTS']
        v_['INCR'] = v_['ICR0']*v_['2PI']
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['THETA1'] = v_['RI-1']*v_['INCR']
            v_['ST1'] = np.sin(v_['THETA1'])
            v_['CT1'] = np.cos(v_['THETA1'])
            v_['P5CT1'] = v_['TP5']*v_['CT1']
            v_['P4P5CT1'] = v_['TP4']+v_['P5CT1']
            v_['R3'] = v_['TP5']*v_['ST1']
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['THETA2'] = v_['RJ-1']*v_['INCR']
                v_['ST2'] = np.sin(v_['THETA2'])
                v_['CT2'] = np.cos(v_['THETA2'])
                v_['R1'] = v_['P4P5CT1']*v_['CT2']
                v_['R2'] = v_['P4P5CT1']*v_['ST2']
                v_['XSEED'] = v_['THETA2']*v_['PSEED']
                v_['SSEED'] = np.cos(v_['XSEED'])
                v_['PER-1'] = v_['PSIZE']*v_['SSEED']
                v_['PERT'] = 1.0+v_['PER-1']
                v_['XD'+str(I)+','+str(J)] = v_['R1']*v_['PERT']
                v_['YD'+str(I)+','+str(J)] = v_['R2']*v_['PERT']
                v_['ZD'+str(I)+','+str(J)] = v_['R3']*v_['PERT']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2x_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2x_ii('Y'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(J))
                [iv,ix_,_] = s2x_ii('Z'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Z'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                [ig,ig_,_] = s2x_ii('OX'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('OY'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Y'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('OZ'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Z'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
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
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['OX'+str(I)+','+str(J)],float(v_['XD'+str(I)+','+str(J)])))
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['OY'+str(I)+','+str(J)],float(v_['YD'+str(I)+','+str(J)])))
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['OZ'+str(I)+','+str(J)],float(v_['ZD'+str(I)+','+str(J)])))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['P4']] = 0.001
        pb.xlower[ix_['P5']] = 0.001
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('P1' in ix_):
            pb.x0[ix_['P1']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1']),float(1.0)))
        if('P2' in ix_):
            pb.x0[ix_['P2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2']),float(0.0)))
        if('P3' in ix_):
            pb.x0[ix_['P3']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P3']),float(1.0)))
        if('P4' in ix_):
            pb.x0[ix_['P4']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P4']),float(1.0)))
        if('P5' in ix_):
            pb.x0[ix_['P5']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P5']),float(0.5)))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                if('X'+str(I)+','+str(J) in ix_):
                    pb.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['XD'+str(I)+','+str(J)])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(J)]),float(v_['XD'+str(I)+','+str(J)])))
                if('Y'+str(I)+','+str(J) in ix_):
                    pb.x0[ix_['Y'+str(I)+','+str(J)]] = float(v_['YD'+str(I)+','+str(J)])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y'+str(I)+','+str(J)]),float(v_['YD'+str(I)+','+str(J)])))
                if('Z'+str(I)+','+str(J) in ix_):
                    pb.x0[ix_['Z'+str(I)+','+str(J)]] = float(v_['ZD'+str(I)+','+str(J)])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Z'+str(I)+','+str(J)]),float(v_['ZD'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eTA', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        elftv = loaset(elftv,it,2,'A')
        elftv = loaset(elftv,it,3,'B')
        elftv = loaset(elftv,it,4,'C')
        [it,iet_,_] = s2x_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        elftv = loaset(elftv,it,1,'P')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'XX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                ename = 'EA'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eTA')
                ielftype = arrset(ielftype, ie, iet_["eTA"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'P1'
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='A')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'P2'
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'P4'
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'EB'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'P3'
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='P')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'EC'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'P5'
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
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
            for J in range(int(v_['1']),int(v_['NPTS'])+1):
                ig = ig_['OX'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL2')
                ig = ig_['OY'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL2')
                ig = ig_['OZ'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL2')
                ig = ig_['A'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EC'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
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
        CCSQ = IV_[2]*IV_[2]
        CCCB = CCSQ*IV_[2]
        XXYY = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        T = XXYY/CCSQ
        DTDX = 2.0*IV_[0]/CCSQ
        DTDY = 2.0*IV_[1]/CCSQ
        DTDC = -2.0*XXYY/CCCB
        D2TDX2 = 2.0/CCSQ
        D2TDY2 = 2.0/CCSQ
        D2TDC2 = 6.0*XXYY/(CCSQ*CCSQ)
        D2TDXC = -4.0*IV_[0]/CCCB
        D2TDYC = -4.0*IV_[1]/CCCB
        S = np.sqrt(T)
        R = 0.5/S
        DSDX = R*DTDX
        DSDY = R*DTDY
        DSDC = R*DTDC
        D2SDX2 = R*(D2TDX2-0.5*DTDX*DTDX/T)
        D2SDY2 = R*(D2TDY2-0.5*DTDY*DTDY/T)
        D2SDC2 = R*(D2TDC2-0.5*DTDC*DTDC/T)
        D2SDXY = -0.5*DTDX*DSDY/T
        D2SDXC = R*(D2TDXC-0.5*DTDX*DTDC/T)
        D2SDYC = R*(D2TDYC-0.5*DTDY*DTDC/T)
        SS = S-1.0
        SPS = SS+SS
        f_   = SS*SS
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SPS*DSDX
            g_[1] = SPS*DSDY
            g_[2] = SPS*DSDC
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = SPS*D2SDX2+2.0*DSDX*DSDX
                H_[0,1] = SPS*D2SDXY+2.0*DSDX*DSDY
                H_[1,0] = H_[0,1]
                H_[0,2] = SPS*D2SDXC+2.0*DSDX*DSDC
                H_[2,0] = H_[0,2]
                H_[1,1] = SPS*D2SDY2+2.0*DSDY*DSDY
                H_[1,2] = SPS*D2SDYC+2.0*DSDY*DSDC
                H_[2,1] = H_[1,2]
                H_[2,2] = SPS*D2SDC2+2.0*DSDC*DSDC
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eISQ(pbm,nargout,*args):

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

