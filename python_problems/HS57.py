from s2mpjlib import *
class  HS57(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS57
#    *********
# 
#    Source: problem 57 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, April 1990
# 
#    classification = "SQR2-AN-2-1"
# 
#    Problem parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS57'

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
        v_['A1'] = 8.0
        v_['A2'] = 8.0
        v_['A3'] = 10.0
        v_['A4'] = 10.0
        v_['A5'] = 10.0
        v_['A6'] = 10.0
        v_['A7'] = 12.0
        v_['A8'] = 12.0
        v_['A9'] = 12.0
        v_['A10'] = 12.0
        v_['A11'] = 14.0
        v_['A12'] = 14.0
        v_['A13'] = 14.0
        v_['A14'] = 16.0
        v_['A15'] = 16.0
        v_['A16'] = 16.0
        v_['A17'] = 18.0
        v_['A18'] = 18.0
        v_['A19'] = 20.0
        v_['A20'] = 20.0
        v_['A21'] = 20.0
        v_['A22'] = 22.0
        v_['A23'] = 22.0
        v_['A24'] = 22.0
        v_['A25'] = 24.0
        v_['A26'] = 24.0
        v_['A27'] = 24.0
        v_['A28'] = 26.0
        v_['A29'] = 26.0
        v_['A30'] = 26.0
        v_['A31'] = 28.0
        v_['A32'] = 28.0
        v_['A33'] = 30.0
        v_['A34'] = 30.0
        v_['A35'] = 30.0
        v_['A36'] = 32.0
        v_['A37'] = 32.0
        v_['A38'] = 34.0
        v_['A39'] = 36.0
        v_['A40'] = 36.0
        v_['A41'] = 38.0
        v_['A42'] = 38.0
        v_['A43'] = 40.0
        v_['A44'] = 42.0
        v_['B1'] = 0.49
        v_['B2'] = 0.49
        v_['B3'] = 0.48
        v_['B4'] = 0.47
        v_['B5'] = 0.48
        v_['B6'] = 0.47
        v_['B7'] = 0.46
        v_['B8'] = 0.46
        v_['B9'] = 0.45
        v_['B10'] = 0.43
        v_['B11'] = 0.45
        v_['B12'] = 0.43
        v_['B13'] = 0.43
        v_['B14'] = 0.44
        v_['B15'] = 0.43
        v_['B16'] = 0.43
        v_['B17'] = 0.46
        v_['B18'] = 0.45
        v_['B19'] = 0.42
        v_['B20'] = 0.42
        v_['B21'] = 0.43
        v_['B22'] = 0.41
        v_['B23'] = 0.41
        v_['B24'] = 0.40
        v_['B25'] = 0.42
        v_['B26'] = 0.40
        v_['B27'] = 0.40
        v_['B28'] = 0.41
        v_['B29'] = 0.40
        v_['B30'] = 0.41
        v_['B31'] = 0.41
        v_['B32'] = 0.40
        v_['B33'] = 0.40
        v_['B34'] = 0.40
        v_['B35'] = 0.38
        v_['B36'] = 0.41
        v_['B37'] = 0.40
        v_['B38'] = 0.40
        v_['B39'] = 0.41
        v_['B40'] = 0.38
        v_['B41'] = 0.40
        v_['B42'] = 0.40
        v_['B43'] = 0.39
        v_['B44'] = 0.39
        v_['1'] = 1
        v_['44'] = 44
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(0.49)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(0.09))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 0.4
        pb.xlower[ix_['X2']] = -4.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(0.42)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(0.42)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(5.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(5.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'AA')
        elftp = loaset(elftp,it,1,'BB')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['44'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eOBSQ')
            ielftype = arrset(ielftype, ie, iet_["eOBSQ"])
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='AA')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(I)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='BB')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(I)]))
        ename = 'PR'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['44'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['CON1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PR'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                 0.02845966
# LO SOLTN                 0.03063791
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
        pb.pbclass = "SQR2-AN-2-1"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOBSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        AM8 = pbm.elpar[iel_][0]-8.0
        CMV1 = 0.49-EV_[0]
        E = np.exp(-EV_[1]*AM8)
        DED2 = -AM8*E
        R = pbm.elpar[iel_][1]-EV_[0]-CMV1*E
        DRD1 = E-1.0
        DRD2 = -CMV1*DED2
        D2RD22 = -CMV1*AM8*AM8*E
        f_   = R*R
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*R*DRD1
            g_[1] = 2.0*R*DRD2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*DRD1*DRD1
                H_[0,1] = 2.0*(DRD2*DRD1+R*DED2)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*(DRD2*DRD2+R*D2RD22)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(pbm,nargout,*args):

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

