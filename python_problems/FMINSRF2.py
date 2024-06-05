from s2mpjlib import *
class  FMINSRF2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FMINSRF2
#    *********
# 
#    The free boundary minimum surface problem.
# 
#    The problem comes from the discretization of the minimum surface
#    problem on the unit square with "free boundary conditions"
#    one must find the minumum surface over the unit square 
#    (which is clearly 1.0).  Furthermore, the distance of the surface
#    from zero at the centre of the unit square is also minimized.
# 
#    The unit square is discretized into (p-1)**2 little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are p**2 of them.
#    Given these heights, the area above a little square is
#    approximated by the
#      S(i,j) = sqrt( 1 + 0.5(p-1)**2 ( a(i,j) + b(i,j) ) ) / (p-1)**2
#    where
#      a(i,j) = x(i,j) - x(i+1,j+1)
#    and
#      b(i,j) = x(i+1,j) - x(i,j+1)
# 
#    Source: setting the boundary free in 
#    A Griewank and Ph. Toint,
#    "Partitioned variable metric updates for large structured
#    optimization problems",
#    Numerische Mathematik 39:429-448, 1982.
# 
#    SIF input: Ph. Toint, November 1991.
# 
#    classification = "OUR2-MY-V-0"
# 
#    P is the number of points in one side of the unit square
# 
#           Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 16
# IE P                   7              $-PARAMETER n = 49
# IE P                   8              $-PARAMETER n = 64     original value
# IE P                   11             $-PARAMETER n = 121
# IE P                   31             $-PARAMETER n = 961
# IE P                   32             $-PARAMETER n = 1024
# IE P                   75             $-PARAMETER n = 5625
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FMINSRF2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'FMINSRF2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(4);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   100            $-PARAMETER n = 10000
# IE P                   125            $-PARAMETER n = 15625
        v_['H00'] = 1.0
        v_['SLOPEJ'] = 4.0
        v_['SLOPEI'] = 8.0
        v_['TWOP'] = v_['P']+v_['P']
        v_['P-1'] = -1+v_['P']
        v_['PP-1'] = v_['P']*v_['P-1']
        v_['RP-1'] = float(v_['P-1'])
        v_['INVP-1'] = 1.0/v_['RP-1']
        v_['RP-1SQ'] = v_['INVP-1']*v_['INVP-1']
        v_['SCALE'] = 1.0/v_['RP-1SQ']
        v_['SQP-1'] = v_['RP-1']*v_['RP-1']
        v_['PARAM'] = 0.5*v_['SQP-1']
        v_['RP'] = float(v_['P'])
        v_['P2'] = v_['RP']*v_['RP']
        v_['P4'] = v_['P2']*v_['P2']
        v_['1'] = 1
        v_['2'] = 2
        v_['STON'] = v_['INVP-1']*v_['SLOPEI']
        v_['WTOE'] = v_['INVP-1']*v_['SLOPEJ']
        v_['H01'] = v_['H00']+v_['SLOPEJ']
        v_['H10'] = v_['H00']+v_['SLOPEI']
        v_['MID'] = int(np.fix(v_['P']/v_['2']))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                pbm.gscale = arrset(pbm.gscale,ig,float(v_['SCALE']))
        [ig,ig_,_] = s2mpj_ii('MID',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['MID']))+','+str(int(v_['MID']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['P2']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),-1.0)
        pbm.gconst = arrset(pbm.gconst,ig_['MID'],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        for J in range(int(v_['1']),int(v_['P'])+1):
            v_['J-1'] = -1+J
            v_['RJ-1'] = float(v_['J-1'])
            v_['TH'] = v_['RJ-1']*v_['WTOE']
            v_['TL'] = v_['TH']+v_['H00']
            v_['TU'] = v_['TH']+v_['H10']
            pb.x0[ix_['X'+str(int(v_['1']))+','+str(J)]] = float(v_['TL'])
            pb.x0[ix_['X'+str(int(v_['P']))+','+str(J)]] = float(v_['TU'])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TV'] = v_['RI-1']*v_['STON']
            v_['TR'] = v_['TV']+v_['H00']
            v_['TL'] = v_['TV']+v_['H01']
            pb.x0[ix_['X'+str(I)+','+str(int(v_['P']))]] = float(v_['TL'])
            pb.x0[ix_['X'+str(I)+','+str(int(v_['1']))]] = float(v_['TR'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQROOT',igt_)
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['MID']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                ig = ig_['S'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gSQROOT')
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['PARAM']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['PARAM']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-MY-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQROOT(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        SQRAL = np.sqrt(GVAR_)
        f_= SQRAL
        if nargout>1:
            g_ = 0.5e0/SQRAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.25e0/(SQRAL*GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

