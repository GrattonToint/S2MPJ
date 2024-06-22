from s2mpjlib import *
class  SCOND1LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SCOND1LS
#    *********
# 
#    The semiconductor problem by Rheinboldt, using a finite difference
#    approximation.
#    This is the least-squares version of problem SEMICON1.
# 
#    Source: problem 10 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SBR2-AN-V-V"
# 
#    N  = Number of discretized point inside the interval [a, b]
#    LN = Index of the last negative discretization point
#         (the interest is in the negative part)
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SCOND1LS'

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
        if nargin<2:
            v_['LN'] = int(9);  #  SIF file default value
        else:
            v_['LN'] = int(args[1])
# IE N                   50             $-PARAMETER
# IE LN                  45             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE LN                  90             $-PARAMETER
# IE N                   500            $-PARAMETER
# IE LN                  450            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE LN                  900            $-PARAMETER
# IE N                   5000           $-PARAMETER
# IE LN                  4500           $-PARAMETER
        if nargin<3:
            v_['LAMBDA'] = float(1.0);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[2])
        v_['A'] = -0.00009
        v_['B'] = 0.00001
        v_['UA'] = 0.0
        v_['UB'] = 700.0
        v_['CA'] = 1.0e12
        v_['CB'] = 1.0e13
        v_['BETA'] = 40.0
        v_['LN+1'] = 1+v_['LN']
        v_['N+1'] = 1+v_['N']
        v_['-A'] = -1.0*v_['A']
        v_['B-A'] = v_['B']+v_['-A']
        v_['RN+1'] = float(v_['N+1'])
        v_['TMP'] = 1.0/v_['RN+1']
        v_['H'] = v_['B-A']*v_['TMP']
        v_['H2'] = v_['H']*v_['H']
        v_['LB'] = v_['LAMBDA']*v_['BETA']
        v_['H2CA'] = v_['H2']*v_['CA']
        v_['H2CB'] = v_['H2']*v_['CB']
        v_['LH2CA'] = v_['LAMBDA']*v_['H2CA']
        v_['LH2CB'] = v_['LAMBDA']*v_['H2CB']
        v_['LUA'] = v_['LAMBDA']*v_['UA']
        v_['LUB'] = v_['LAMBDA']*v_['UB']
        v_['ULW'] = -5.0+v_['LUA']
        v_['UUP'] = 5.0+v_['LUB']
        v_['-LB'] = -1.0*v_['LB']
        v_['-LUB'] = -1.0*v_['LUB']
        v_['-LH2CB'] = -1.0*v_['LH2CB']
        v_['0'] = 0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N+1'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['LN'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['LH2CA']))
        for I in range(int(v_['LN+1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['-LH2CB']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),v_['UUP'])
        pb.xlower = np.full((pb.n,1),v_['ULW'])
        pb.xlower[ix_['U'+str(int(v_['0']))]] = v_['LUA']
        pb.xupper[ix_['U'+str(int(v_['0']))]] = v_['LUA']
        pb.xlower[ix_['U'+str(int(v_['N+1']))]] = v_['LUB']
        pb.xupper[ix_['U'+str(int(v_['N+1']))]] = v_['LUB']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        pb.x0[ix_['U'+str(int(v_['0']))]] = float(v_['LUA'])
        pb.x0[ix_['U'+str(int(v_['N+1']))]] = float(v_['LUB'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eWE1', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'LAC')
        elftp = loaset(elftp,it,1,'LAB')
        elftp = loaset(elftp,it,2,'LU')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWE1')
            ielftype = arrset(ielftype, ie, iet_["eWE1"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,v_['ULW'],v_['UUP'],0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='LAC')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['LH2CA']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='LAB')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-LB']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='LU')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['LUA']))
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eWE1')
            ielftype = arrset(ielftype, ie, iet_["eWE1"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,v_['ULW'],v_['UUP'],0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='LAC')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-LH2CB']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='LAB')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['LB']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='LU')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['LUB']))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eWE1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVAL  = (
              pbm.elpar[iel_][0]*np.exp(pbm.elpar[iel_][1]*(EV_[0]-pbm.elpar[iel_][2])))
        f_   = FVAL
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][1]*FVAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*FVAL
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

