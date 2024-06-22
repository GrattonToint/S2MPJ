from s2mpjlib import *
class  HS105(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS105
#    *********
# 
#    Source: problem 105 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
#    bug correction (line 351) Ph. Toint, May 2024
# 
#    classification = "OLR2-AY-8-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS105'

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
        v_['N'] = 8
        v_['1'] = 1
        v_['235'] = 235
        v_['Y1'] = 95.0
        v_['Y2'] = 105.0
        v_['LOW'] = 3
        v_['UP'] = 6
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 110.0
        v_['LOW'] = 7
        v_['UP'] = 10
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 115.0
        v_['LOW'] = 11
        v_['UP'] = 25
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 120.0
        v_['LOW'] = 26
        v_['UP'] = 40
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 125.0
        v_['LOW'] = 41
        v_['UP'] = 55
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 130.0
        v_['LOW'] = 56
        v_['UP'] = 68
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 135.0
        v_['LOW'] = 69
        v_['UP'] = 89
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 140.0
        v_['LOW'] = 90
        v_['UP'] = 101
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 145.0
        v_['LOW'] = 102
        v_['UP'] = 118
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 150.0
        v_['LOW'] = 119
        v_['UP'] = 122
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 155.0
        v_['LOW'] = 123
        v_['UP'] = 142
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 160.0
        v_['LOW'] = 143
        v_['UP'] = 150
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 165.0
        v_['LOW'] = 151
        v_['UP'] = 167
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 170.0
        v_['LOW'] = 168
        v_['UP'] = 175
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 175.0
        v_['LOW'] = 176
        v_['UP'] = 181
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 180.0
        v_['LOW'] = 182
        v_['UP'] = 187
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 185.0
        v_['LOW'] = 188
        v_['UP'] = 194
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 190.0
        v_['LOW'] = 195
        v_['UP'] = 198
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 195.0
        v_['LOW'] = 199
        v_['UP'] = 201
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 200.0
        v_['LOW'] = 202
        v_['UP'] = 204
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 205.0
        v_['LOW'] = 205
        v_['UP'] = 212
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 210.0
        v_['Y213'] = 215.0
        v_['LOW'] = 214
        v_['UP'] = 219
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 220.0
        v_['LOW'] = 220
        v_['UP'] = 224
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 230.0
        v_['Y225'] = 235.0
        v_['LOW'] = 226
        v_['UP'] = 232
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 240.0
        v_['Y233'] = 245.0
        v_['LOW'] = 234
        v_['UP'] = 235
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 250.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['235'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0e+0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['C1'],float(-1.0e+0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 0.001
        pb.xupper[ix_['X1']] = 0.499
        pb.xlower[ix_['X2']] = 0.001
        pb.xupper[ix_['X2']] = 0.499
        pb.xlower[ix_['X3']] = 100.0
        pb.xupper[ix_['X3']] = 180.0
        pb.xlower[ix_['X4']] = 130.0
        pb.xupper[ix_['X4']] = 210.0
        pb.xlower[ix_['X5']] = 170.0
        pb.xupper[ix_['X5']] = 240.0
        pb.xlower[ix_['X6']] = 5.0
        pb.xupper[ix_['X6']] = 25.0
        pb.xlower[ix_['X7']] = 5.0
        pb.xupper[ix_['X7']] = 25.0
        pb.xlower[ix_['X8']] = 5.0
        pb.xupper[ix_['X8']] = 25.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(0.1)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(0.2)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),float(0.2))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(100.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(100.0)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(125.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),float(125.0))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(175.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(175.0)))
        if('X6' in ix_):
            pb.x0[ix_['X6']] = float(11.2)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X6']),float(11.2))
        if('X7' in ix_):
            pb.x0[ix_['X7']] = float(13.2)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X7']),float(13.2)))
        if('X8' in ix_):
            pb.x0[ix_['X8']] = float(15.8)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X8']),float(15.8))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eABI', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'YI')
        [it,iet_,_] = s2mpj_ii( 'eCI', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X5')
        elftv = loaset(elftv,it,3,'X8')
        elftp = loaset(elftp,it,0,'YI')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['235'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eABI')
            ielftype = arrset(ielftype, ie, iet_["eABI"])
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X6'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='YI')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eABI')
            ielftype = arrset(ielftype, ie, iet_["eABI"])
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X7'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='YI')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCI')
            ielftype = arrset(ielftype, ie, iet_["eCI"])
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X8'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X8')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='YI')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['235'])+1):
            ig = ig_['OBJ'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-AY-8-1"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eABI(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R = EV_[0]/EV_[1]
        D = (pbm.elpar[iel_][0]-EV_[2])/EV_[1]
        E = np.exp(-5.0e-1*D*D)
        DDV2 = -D/EV_[1]
        DEV2 = E*(-D)*DDV2
        DDV3 = -1.0e+0/EV_[1]
        DEV3 = E*(-D)*DDV3
        f_   = R*E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/EV_[1]
            g_[1] = (D*D-1.0e+0)*R*E/EV_[1]
            g_[2] = D*R*E/EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEV2-E/EV_[1])/EV_[1]
                H_[1,0] = H_[0,1]
                H_[0,2] = DEV3/EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/EV_[1]))*R/EV_[1]
                H_[1,2] = (DDV2*E+D*DEV2-2.0e+0*D*E/EV_[1])*R/EV_[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = (DDV3*E+D*DEV3)*R/EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCI(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]-1
        U_[0,1] = U_[0,1]-1
        U_[1,3] = U_[1,3]+1
        U_[2,2] = U_[2,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        R = (1.0e+0+IV_[0])/IV_[1]
        D = (pbm.elpar[iel_][0]-IV_[2])/IV_[1]
        E = np.exp(-5.0e-1*D*D)
        DDV2 = -D/IV_[1]
        DEV2 = E*(-D)*DDV2
        DDV3 = -1.0e+0/IV_[1]
        DEV3 = E*(-D)*DDV3
        f_   = R*E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/IV_[1]
            g_[1] = (D*D-1.0e+0)*R*E/IV_[1]
            g_[2] = D*R*E/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEV2-E/IV_[1])/IV_[1]
                H_[1,0] = H_[0,1]
                H_[0,2] = DEV3/IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/IV_[1]))*R/IV_[1]
                H_[1,2] = (DDV2*E+D*DEV2-2.0e+0*D*E/IV_[1])*R/IV_[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = (DDV3*E+D*DEV3)*R/IV_[1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def g_globs(pbm):

        pbm.gfpar = np.array([]);
        pbm.gfpar = arrset( pbm.gfpar,0,3.9894228040143270e-01)    # this is P
        return pbm

    @staticmethod
    def gLOG(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= -np.log(pbm.gfpar[0]*GVAR_)
        if nargout>1:
            g_ = -1/GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 1/GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

