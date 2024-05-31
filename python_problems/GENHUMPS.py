from s2xlib import *
class  GENHUMPS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GENHUMPS
#    *********
# 
#    A multi-dimensional variant of HUMPS, a two dimensional function
#    with a lot of humps. The density of humps increases with the
#    parameter ZETA, making the problem more difficult.
# 
#    The problem is nonconvex.
# 
#    Source:
#    Ph. Toint, private communication, 1997.
# 
#    SDIF input: N. Gould and Ph. Toint, November 1997.
# 
#    classification = "OUR2-AN-V-0"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   5              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GENHUMPS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'GENHUMPS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER    original value
# IE N                   5000           $-PARAMETER
# RE ZETA                2.0            $-PARAMETER
        if nargin<2:
            v_['ZETA'] = float(20.0);  #  SIF file default value
        else:
            v_['ZETA'] = float(args[1])
        v_['1'] = 1
        v_['N-1'] = -1+v_['N']
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
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(-506.2))
        pb.x0[ix_['X'+str(int(v_['1']))]] = float(-506.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eHMP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ename = 'Q'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eHMP')
            ielftype = arrset(ielftype, ie, iet_["eHMP"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-506.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-506.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ZETA']))
            ename = 'SX'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-506.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'SY'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-506.2)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SX'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.05))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SY'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.05))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHMP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SAX = np.sin(pbm.elpar[iel_][0]*EV_[0])
        SAY = np.sin(pbm.elpar[iel_][0]*EV_[1])
        CAX = np.cos(pbm.elpar[iel_][0]*EV_[0])
        CAY = np.cos(pbm.elpar[iel_][0]*EV_[1])
        AA = pbm.elpar[iel_][0]+pbm.elpar[iel_][0]
        AAA = pbm.elpar[iel_][0]*AA
        f_   = (SAX*SAY)**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = AA*SAX*CAX*SAY**2
            g_[1] = AA*SAX**2*CAY*SAY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = AAA*SAY**2*(CAX**2-SAX**2)
                H_[0,1] = AA*AA*SAX*CAX*SAY*CAY
                H_[1,0] = H_[0,1]
                H_[1,1] = AAA*SAX**2*(CAY**2-SAY**2)
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
        f_   = EV_[0]**2
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

