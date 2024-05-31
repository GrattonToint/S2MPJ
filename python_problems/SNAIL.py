from s2xlib import *
class  SNAIL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SNAIL
#    *********
# 
#    A 2D problem featuring a spiraling valley.
#    Dedicated to the city of Namur, whose emblem is a snail.
# 
#    Source:
#    J. Engels, private communication.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "OUR2-AN-2-0"
# 
#    Problem parameters (CUP > CLOW > 0)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SNAIL'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'SNAIL'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['CLOW'] = 1.0
        v_['CUP'] = 2.0
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
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(10.0)
        pb.x0[ix_['X2']] = float(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSPIRAL', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'CL')
        elftp = loaset(elftp,it,1,'CU')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSPIRAL')
        ielftype = arrset(ielftype, ie, iet_["eSPIRAL"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='CL')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CLOW']))
        posep = find(elftp[ielftype[ie]],lambda x:x=='CU')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CUP']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-2-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSPIRAL(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        A = 0.5*(pbm.elpar[iel_][1]+pbm.elpar[iel_][0])
        B = 0.5*(pbm.elpar[iel_][1]-pbm.elpar[iel_][0])
        X2 = EV_[0]*EV_[0]
        Y2 = EV_[1]*EV_[1]
        R2 = X2+Y2
        D = 1.0+R2
        D2 = D*D
        D3 = D2*D
        U = R2/D
        DUDX = (EV_[0]+EV_[0])/D2
        DUDY = (EV_[1]+EV_[1])/D2
        D2UDX2 = 2.0*(D-4.0*X2)/D3
        D2UDY2 = 2.0*(D-4.0*Y2)/D3
        D2UDXY = -8.0*EV_[0]*EV_[1]/D3
        THETA = np.arctan2(EV_[1],EV_[0])
        DTDX = -EV_[1]/R2
        DTDY = EV_[0]/R2
        R4 = R2*R2
        D2TDX2 = 2.0*EV_[0]*EV_[1]/R4
        D2TDY2 = -2.0*EV_[1]*EV_[0]/R4
        D2TDXY = (Y2-X2)/R4
        R = np.sqrt(R2)
        R3 = R*R2
        DRDX = EV_[0]/R
        DRDY = EV_[1]/R
        D2RDX2 = Y2/R3
        D2RDY2 = X2/R3
        D2RDXY = -EV_[0]*EV_[1]/R3
        ARG = R-THETA
        S = B*np.sin(ARG)
        C = B*np.cos(ARG)
        DCDX = -S*(DRDX-DTDX)
        DCDY = -S*(DRDY-DTDY)
        D2CDX2 = -C*(DRDX-DTDX)**2-S*(D2RDX2-D2TDX2)
        D2CDY2 = -C*(DRDY-DTDY)**2-S*(D2RDY2-D2TDY2)
        D2CDXY = -C*(DRDX-DTDX)*(DRDY-DTDY)-S*(D2RDXY-D2TDXY)
        V = 1.0+A*R-R*C
        DVDX = A*DRDX-DRDX*C-R*DCDX
        DVDY = A*DRDY-DRDY*C-R*DCDY
        D2VDX2 = A*D2RDX2-D2RDX2*C-2.0*DRDX*DCDX-R*D2CDX2
        D2VDY2 = A*D2RDY2-D2RDY2*C-2.0*DRDY*DCDY-R*D2CDY2
        D2VDXY = A*D2RDXY-D2RDXY*C-DRDX*DCDY-DRDY*DCDX-R*D2CDXY
        f_   = U*V
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = DUDX*V+U*DVDX
            g_[1] = DUDY*V+U*DVDY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = D2UDX2*V+2.0*DUDX*DVDX+U*D2VDX2
                H_[0,1] = D2UDXY*V+DUDX*DVDY+DUDY*DVDX+U*D2VDXY
                H_[1,0] = H_[0,1]
                H_[1,1] = D2UDY2*V+2.0*DUDY*DVDY+U*D2VDY2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

