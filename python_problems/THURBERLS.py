from s2xlib import *
class  THURBERLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : THURBERLS
#    *********
# 
#    NIST Data fitting problem THURBERLS.
# 
#    Fit: y = (b1 + b2*x + b3*x**2 + b4*x**3) / 
#             (1 + b5*x + b6*x**2 + b7*x**3) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Thurber, R., NIST (197?).  
#      Semiconductor electron mobility modeling.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "SUR2-MN-7-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'THURBERLS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'THURBERLS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 37
        v_['N'] = 7
        v_['1'] = 1
        v_['X1'] = -3.067
        v_['X2'] = -2.981
        v_['X3'] = -2.921
        v_['X4'] = -2.912
        v_['X5'] = -2.840
        v_['X6'] = -2.797
        v_['X7'] = -2.702
        v_['X8'] = -2.699
        v_['X9'] = -2.633
        v_['X10'] = -2.481
        v_['X11'] = -2.363
        v_['X12'] = -2.322
        v_['X13'] = -1.501
        v_['X14'] = -1.460
        v_['X15'] = -1.274
        v_['X16'] = -1.212
        v_['X17'] = -1.100
        v_['X18'] = -1.046
        v_['X19'] = -0.915
        v_['X20'] = -0.714
        v_['X21'] = -0.566
        v_['X22'] = -0.545
        v_['X23'] = -0.400
        v_['X24'] = -0.309
        v_['X25'] = -0.109
        v_['X26'] = -0.103
        v_['X27'] = 0.010
        v_['X28'] = 0.119
        v_['X29'] = 0.377
        v_['X30'] = 0.790
        v_['X31'] = 0.963
        v_['X32'] = 1.006
        v_['X33'] = 1.115
        v_['X34'] = 1.572
        v_['X35'] = 1.841
        v_['X36'] = 2.047
        v_['X37'] = 2.200
        v_['Y1'] = 80.574
        v_['Y2'] = 84.248
        v_['Y3'] = 87.264
        v_['Y4'] = 87.195
        v_['Y5'] = 89.076
        v_['Y6'] = 89.608
        v_['Y7'] = 89.868
        v_['Y8'] = 90.101
        v_['Y9'] = 92.405
        v_['Y10'] = 95.854
        v_['Y11'] = 100.696
        v_['Y12'] = 101.060
        v_['Y13'] = 401.672
        v_['Y14'] = 390.724
        v_['Y15'] = 567.534
        v_['Y16'] = 635.316
        v_['Y17'] = 733.054
        v_['Y18'] = 759.087
        v_['Y19'] = 894.206
        v_['Y20'] = 990.785
        v_['Y21'] = 1090.109
        v_['Y22'] = 1080.914
        v_['Y23'] = 1122.643
        v_['Y24'] = 1178.351
        v_['Y25'] = 1260.531
        v_['Y26'] = 1273.514
        v_['Y27'] = 1288.339
        v_['Y28'] = 1327.543
        v_['Y29'] = 1353.863
        v_['Y30'] = 1414.509
        v_['Y31'] = 1425.208
        v_['Y32'] = 1421.384
        v_['Y33'] = 1442.962
        v_['Y34'] = 1464.350
        v_['Y35'] = 1468.705
        v_['Y36'] = 1447.894
        v_['Y37'] = 1457.628
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2x_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['B1']] = float(1000.0)
        pb.x0[ix_['B2']] = float(1000.0)
        pb.x0[ix_['B3']] = float(400.0)
        pb.x0[ix_['B4']] = float(40.0)
        pb.x0[ix_['B5']] = float(0.7)
        pb.x0[ix_['B6']] = float(0.3)
        pb.x0[ix_['B7']] = float(0.03)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE19', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE19')
            ielftype = arrset(ielftype, ie, iet_["eE19"])
            vname = 'B1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B5'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B6'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B7'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V7')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['X'+str(I)]))
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
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-MN-7-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE19(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        X2 = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]
        X3 = X2*pbm.elpar[iel_][0]
        X4 = X3*pbm.elpar[iel_][0]
        X5 = X4*pbm.elpar[iel_][0]
        X6 = X5*pbm.elpar[iel_][0]
        T = EV_[0]+EV_[1]*pbm.elpar[iel_][0]+EV_[2]*X2+EV_[3]*X3
        D = 1.0e0+EV_[4]*pbm.elpar[iel_][0]+EV_[5]*X2+EV_[6]*X3
        D2 = D*D
        TD3 = 0.5e0*D2*D
        f_   = T/D
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0e0/D
            g_[1] = pbm.elpar[iel_][0]/D
            g_[2] = X2/D
            g_[3] = X3/D
            g_[4] = -pbm.elpar[iel_][0]*T/D2
            g_[5] = -X2*T/D2
            g_[6] = -X3*T/D2
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,4] = -pbm.elpar[iel_][0]/D2
                H_[4,0] = H_[0,4]
                H_[0,5] = -X2/D2
                H_[5,0] = H_[0,5]
                H_[0,6] = -X3/D2
                H_[6,0] = H_[0,6]
                H_[1,4] = -X2/D2
                H_[4,1] = H_[1,4]
                H_[1,5] = -X3/D2
                H_[5,1] = H_[1,5]
                H_[1,6] = -X4/D2
                H_[6,1] = H_[1,6]
                H_[2,4] = -X3/D2
                H_[4,2] = H_[2,4]
                H_[2,5] = -X4/D2
                H_[5,2] = H_[2,5]
                H_[2,6] = -X5/D2
                H_[6,2] = H_[2,6]
                H_[3,4] = -X4/D2
                H_[4,3] = H_[3,4]
                H_[3,5] = -X5/D2
                H_[5,3] = H_[3,5]
                H_[3,6] = -X6/D2
                H_[6,3] = H_[3,6]
                H_[4,4] = X2*T/TD3
                H_[4,5] = X3*T/TD3
                H_[5,4] = H_[4,5]
                H_[4,6] = X4*T/TD3
                H_[6,4] = H_[4,6]
                H_[5,5] = X4*T/TD3
                H_[5,6] = X5*T/TD3
                H_[6,5] = H_[5,6]
                H_[6,6] = X6*T/TD3
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

