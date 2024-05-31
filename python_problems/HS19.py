from s2xlib import *
class  HS19(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS19
#    *********
# 
#    Source: problem 19 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn April 1990
# 
#    classification = "OQR2-AN-2-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS19'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS19'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2x_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
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
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        [ig,ig_,_] = s2x_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON2')
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
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(100.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON2'],float(-82.81))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 13.0
        pb.xupper[ix_['X1']] = 100.0
        pb.xupper[ix_['X2']] = 100.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(20.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(20.1)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(5.84)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(5.84)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eCBm10', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2x_ii( 'eCBm20', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2x_ii( 'eSQm5', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2x_ii( 'eSQm6', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCBm10')
        ielftype = arrset(ielftype, ie, iet_["eCBm10"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCBm20')
        ielftype = arrset(ielftype, ie, iet_["eCBm20"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype, ie, iet_["eSQm5"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm5')
        ielftype = arrset(ielftype, ie, iet_["eSQm5"])
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQm6')
        ielftype = arrset(ielftype, ie, iet_["eSQm6"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
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
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['CON1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['CON2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OQR2-AN-2-2"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCBm10(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-10.0
        f_   = DIF**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCBm20(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-20.0
        f_   = DIF**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*DIF*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*DIF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQm5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-5.0
        f_   = DIF**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQm6(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-6.0
        f_   = DIF**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
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

