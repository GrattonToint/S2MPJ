from s2mpjlib import *
class  POWELLBC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POWELLBC
#    --------
# 
#    A bound-constrained optimization problem to 
#    separate points within a square in the plane
# 
#    Given points p_j = ( x_2j-1 , x_2j ), i = 1, ..., p
# 
#    minimize sum_k=2^p sum_j=1^k-1 1 / || p_j - p_k ||_2
# 
#    subject to 0 <= x_i <= 1, i = 1, ..., 2p = n
# 
#    Source: 
#    M. J. D. Powell
#    Private communication (Optbridge, 2006)
# 
#    SIF input: Nick Gould, Aug 2006.
# 
#    classification = "OBR2-AN-V-0"
# 
#    Number of points
# 
#           Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER
# IE P                   5              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'POWELLBC'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'POWELLBC'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(12);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   10             $-PARAMETER
# IE P                   100            $-PARAMETER
# IE P                   500            $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['N'] = 2*v_['P']
        v_['RN'] = float(v_['N'])
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
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),1.0)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'] = v_['RI']/v_['RN']
            v_['T'] = v_['T']*v_['T']
            pb.x0[ix_['X'+str(I)]] = float(v_['T'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eINVNRM', iet_)
        elftv = loaset(elftv,it,0,'XJ')
        elftv = loaset(elftv,it,1,'YJ')
        elftv = loaset(elftv,it,2,'XK')
        elftv = loaset(elftv,it,3,'YK')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for K in range(int(v_['2']),int(v_['P'])+1):
            v_['K-1'] = K-v_['1']
            v_['2K'] = v_['2']*K
            v_['2K-1'] = v_['2K']-v_['1']
            for J in range(int(v_['1']),int(v_['K-1'])+1):
                v_['2J'] = v_['2']*J
                v_['2J-1'] = v_['2J']-v_['1']
                ename = 'E'+str(K)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eINVNRM')
                    ielftype = arrset( ielftype,ie,iet_['eINVNRM'])
                vname = 'X'+str(int(v_['2J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XJ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['2K-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XK')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['2J']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YJ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['2K']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.0,1.0,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YK')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['2']),int(v_['P'])+1):
            v_['K-1'] = K-v_['1']
            for J in range(int(v_['1']),int(v_['K-1'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(K)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               ??
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-AN-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eINVNRM(pbm,nargout,*args):

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
        NORM = 1.0/np.sqrt(IV_[0]*IV_[0]+IV_[1]*IV_[1])
        f_   = NORM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -IV_[0]*NORM**3
            g_[1] = -IV_[1]*NORM**3
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = (3.0*IV_[0]*IV_[0]*NORM**2-1.0)*NORM**3
                H_[0,1] = 3.0*IV_[0]*IV_[1]*NORM**5
                H_[1,0] = H_[0,1]
                H_[1,1] = (3.0*IV_[1]*IV_[1]*NORM**2-1.0)*NORM**3
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

