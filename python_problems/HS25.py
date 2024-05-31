from s2xlib import *
class  HS25(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS25
#    *********
# 
#    A nonlinear least squares problem with bounds.
# 
#    Source: problem 25 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: J-M Collin, Mar 1990.
# 
#    classification = "SBR2-AN-3-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS25'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS25'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        v_['1'] = 1
        v_['99'] = 99
        v_['2/3'] = 0.6666666666
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2x_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2x_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['99'])+1):
            [ig,ig_,_] = s2x_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['99'])+1):
            v_['IR'] = float(I)
            v_['I/100'] = 0.01*v_['IR']
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)],float(v_['I/100']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X1']] = 0.1
        pb.xupper[ix_['X1']] = 100.0
        pb.xlower[ix_['X2']] = 0.0
        pb.xupper[ix_['X2']] = 25.6
        pb.xlower[ix_['X3']] = 0.0
        pb.xupper[ix_['X3']] = 5.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(100.0)
        pb.x0[ix_['X2']] = float(12.5)
        pb.x0[ix_['X3']] = float(3.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eWFI', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        elftp = []
        elftp = loaset(elftp,it,0,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['99'])+1):
            v_['IR'] = float(I)
            v_['I/100'] = 0.01*v_['IR']
            v_['LOG01I'] = np.log(v_['I/100'])
            v_['M50LOG'] = -50.0*v_['LOG01I']
            v_['EXPLOG'] = np.log(v_['M50LOG'])
            v_['EXPL2/3'] = v_['EXPLOG']*v_['2/3']
            v_['EXP2/3'] = np.exp(v_['EXPL2/3'])
            v_['UI'] = 25.0+v_['EXP2/3']
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWFI')
            ielftype = arrset(ielftype, ie, iet_["eWFI"])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='W')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['UI']))
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
        for I in range(int(v_['1']),int(v_['99'])+1):
            ig = ig_['O'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-AN-3-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eWFI(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XI = 1.0/EV_[0]
        X2I = XI*XI
        X3I = X2I*XI
        WMY = pbm.elpar[iel_][0]-EV_[1]
        WMYEZ = WMY**EV_[2]
        LWMY = np.log(WMY)
        EXPO = np.exp(-XI*WMYEZ)
        f_   = EXPO
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = X2I*WMYEZ*EXPO
            g_[1] = XI*EV_[2]*WMY**(EV_[2]-1.0)*EXPO
            g_[2] = -XI*LWMY*WMYEZ*EXPO
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EXPO*WMYEZ*X3I*(-2.0+XI*WMY**EV_[2])
                H_[0,1] = EXPO*EV_[2]*X2I*WMY**(EV_[2]-1.0)*(-1.0+XI*WMYEZ)
                H_[1,0] = H_[0,1]
                H_[0,2] = EXPO*X2I*WMYEZ*LWMY*(1.0-XI*WMYEZ)
                H_[2,0] = H_[0,2]
                H_[1,1] = EXPO*XI*WMY**(EV_[2]-2.0)*EV_[2]*(-EV_[2]+1.0+XI*EV_[2]*WMYEZ)
                H_[1,2] = EXPO*XI*WMY**(EV_[2]-1.0)*(1.0+EV_[2]*LWMY*(1.0-XI*WMYEZ))
                H_[2,1] = H_[1,2]
                H_[2,2] = EXPO*WMYEZ*XI*LWMY**2*(-1.0+XI*WMYEZ)
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
