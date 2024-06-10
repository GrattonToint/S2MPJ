from s2mpjlib import *
class  S368(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : S368
#    *********
# 
#    Wolfe's problem.
# 
#    Source:
#    P. Wolfe,
#    "Explicit solution of an optimization problem",
#    Mathematical Programming 2, 258-260, 1972.
# 
#    SIF input: Nick Gould, Oct 1992.
# 
#    See also Schittkowski #368 (for N = 8)
# 
#    classification = "OBR2-MN-V-0"
# 
#    The number of variables is N.
# 
#           Alternative values for the SIF file parameters:
# IE N                   8              $-PARAMETER Schittkowski #368
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'S368'

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
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
        v_['1'] = 1
        v_['N+1'] = 1+v_['N']
        v_['RN+1'] = float(v_['N+1'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
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
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('M'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                [ig,ig_,_] = s2mpj_ii('P'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),1.0)
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'] = v_['RI']/v_['RN+1']
            pb.x0[ix_['X'+str(I)]] = float(v_['T'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePRODM', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'ePRODP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ename = 'M'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePRODM')
                ielftype = arrset(ielftype, ie, iet_["ePRODM"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePRODP')
                ielftype = arrset(ielftype, ie, iet_["ePRODP"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['M'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['M'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['P'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               ??
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-MN-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePRODM(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]**2*EV_[1]**4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -2.0e+0*EV_[0]*EV_[1]**4
            g_[1] = -4.0e+0*EV_[0]**2*EV_[1]**3
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = -2.0e+0*EV_[1]**4
                H_[1,0] = -8.0e+0*EV_[0]*EV_[1]**3
                H_[0,1] = H_[1,0]
                H_[1,1] = -1.2e+1*EV_[0]**2*EV_[1]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**3*EV_[1]**3
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*EV_[0]**2*EV_[1]**3
            g_[1] = 3.0e+0*EV_[0]**3*EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 6.0e+0*EV_[0]*EV_[1]**3
                H_[1,0] = 9.0e+0*EV_[0]**2*EV_[1]**2
                H_[0,1] = H_[1,0]
                H_[1,1] = 6.0e+0*EV_[0]**3*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

