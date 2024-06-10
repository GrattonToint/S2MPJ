from s2mpjlib import *
class  DRCAV2LQ(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRCAV2LQ
#    *********
# 
#    This system of nonlinear equations models the stream function
#    corresponding to an incompressible fluid flow in a driven cavity 
#    (after elimination of the vorticity). The system is solved in the
#    least-squares sense.  The nonlinear system formulation is problem DRCAVTY2.
# 
#    The problem is nonconvex.
#    It differs from the problems DRCAV1LQ and DRCAV3LQ by the value 
#    chosen for the Reynolds number.
# 
#    Source:  
#    P.N. Brown and Y. Saad, 
#    "Hybrid Krylov Methods for Nonlinear Systems of Equations",
#    SIAM J. Sci. Stat. Comput. 11, pp. 450-481, 1990.
#    The boundary conditions have been set according to
#    I.E. Kaporin and O. Axelsson,
#    "On a class of nonlinear equation solvers based on the residual norm
#    reduction over a sequence of affine subspaces",
#    SIAM J, Sci. Comput. 16(1), 1995.
# 
#    SIF input: Ph. Toint, Jan 1995.
# 
#    classification = "OXR2-MY-V-V"
# 
#    Discretization mesh: n = (M+3)**2 - fixed variables
# 
#           Alternative values for the SIF file parameters:
# IE M                   10             $-PARAMETER  n =   100   original value
# IE M                   31             $-PARAMETER  n =   961
# IE M                   63             $-PARAMETER  n =  3969
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DRCAV2LQ'

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
            v_['M'] = int(10);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
# IE M                   100            $-PARAMETER  n = 10000
        if nargin<2:
            v_['RE'] = float(1000.0);  #  SIF file default value
        else:
            v_['RE'] = float(args[1])
        v_['M+2'] = 2+v_['M']
        v_['RM+2'] = float(v_['M+2'])
        v_['H'] = 1.0/v_['RM+2']
        v_['-1'] = -1
        v_['0'] = 0
        v_['1'] = 1
        v_['M+1'] = 1+v_['M']
        v_['H/2'] = 0.5*v_['H']
        v_['-H/2'] = -1.0*v_['H/2']
        v_['RE/4'] = 0.25*v_['RE']
        v_['-RE/4'] = -1.0*v_['RE/4']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['-1']),int(v_['M+2'])+1):
            for J in range(int(v_['-1']),int(v_['M+2'])+1):
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['I-2'] = -2+I
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            for J in range(int(v_['1']),int(v_['M'])+1):
                v_['J-2'] = -2+J
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                v_['J+2'] = 2+J
                [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Y'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I-1']))+','+str(J)]
                pbm.A[ig,iv] = float(-8.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+1']))+','+str(J)]
                pbm.A[ig,iv] = float(-8.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(-8.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(-8.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I-1']))+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+1']))+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I-1']))+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+1']))+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I-2']))+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(int(v_['I+2']))+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)+','+str(int(v_['J-2']))]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)+','+str(int(v_['J+2']))]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for J in range(int(v_['-1']),int(v_['M+2'])+1):
            pb.xlower[ix_['Y'+str(int(v_['-1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['Y'+str(int(v_['-1']))+','+str(J)]] = 0.0
            pb.xlower[ix_['Y'+str(int(v_['0']))+','+str(J)]] = 0.0
            pb.xupper[ix_['Y'+str(int(v_['0']))+','+str(J)]] = 0.0
        for I in range(int(v_['1']),int(v_['M'])+1):
            pb.xlower[ix_['Y'+str(I)+','+str(int(v_['-1']))]] = 0.0
            pb.xupper[ix_['Y'+str(I)+','+str(int(v_['-1']))]] = 0.0
            pb.xlower[ix_['Y'+str(I)+','+str(int(v_['0']))]] = 0.0
            pb.xupper[ix_['Y'+str(I)+','+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['M'])+1):
            pb.xlower[ix_['Y'+str(I)+','+str(int(v_['M+1']))]] = 0.0
            pb.xupper[ix_['Y'+str(I)+','+str(int(v_['M+1']))]] = 0.0
            pb.xlower[ix_['Y'+str(I)+','+str(int(v_['M+2']))]] = 0.0
            pb.xupper[ix_['Y'+str(I)+','+str(int(v_['M+2']))]] = 0.0
        for J in range(int(v_['-1']),int(v_['M+2'])+1):
            pb.xlower[ix_['Y'+str(int(v_['M+1']))+','+str(J)]] = v_['-H/2']
            pb.xupper[ix_['Y'+str(int(v_['M+1']))+','+str(J)]] = v_['-H/2']
            pb.xlower[ix_['Y'+str(int(v_['M+2']))+','+str(J)]] = v_['H/2']
            pb.xupper[ix_['Y'+str(int(v_['M+2']))+','+str(J)]] = v_['H/2']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eIPR', iet_)
        elftv = loaset(elftv,it,0,'A1')
        elftv = loaset(elftv,it,1,'A2')
        elftv = loaset(elftv,it,2,'B1')
        elftv = loaset(elftv,it,3,'B2')
        elftv = loaset(elftv,it,4,'B3')
        elftv = loaset(elftv,it,5,'B4')
        elftv = loaset(elftv,it,6,'B5')
        elftv = loaset(elftv,it,7,'B6')
        elftv = loaset(elftv,it,8,'B7')
        elftv = loaset(elftv,it,9,'B8')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['I-2'] = -2+I
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            for J in range(int(v_['1']),int(v_['M'])+1):
                v_['J-2'] = -2+J
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                v_['J+2'] = 2+J
                ename = 'X'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eIPR')
                    ielftype = arrset( ielftype,ie,iet_['eIPR'])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'Y'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-2']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B5')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B6')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B7')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+2']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B8')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Z'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eIPR')
                    ielftype = arrset( ielftype,ie,iet_['eIPR'])
                vname = 'Y'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='A1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='A2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J-2']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B5')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B6')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B7')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J+2']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='B8')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
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
            for J in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['E'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RE/4']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Z'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RE/4']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OXR2-MY-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eIPR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,10))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]-4
        U_[1,6] = U_[1,6]+4
        U_[1,7] = U_[1,7]-1
        U_[1,8] = U_[1,8]-1
        U_[1,9] = U_[1,9]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
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

