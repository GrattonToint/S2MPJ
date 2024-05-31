from s2xlib import *
class  TORSIOND(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TORSIOND
#    *********
# 
#    The quadratic elastic torsion problem
# 
#    The problem comes from the obstacle problem on a square.
# 
#    The square is discretized into (px-1)(py-1) little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are px**2 of them.
# 
#    The dimension of the problem is specified by Q, which is half the
#    number discretization points along one of the coordinate
#    direction.
#    Since the number of variables is P**2, it is given by 4Q**2
# 
#    This is a variant of the problem stated in the report quoted below.
#    It corresponds to the problem as distributed in MINPACK-2.
# 
#    Source: problem (c=10, starting point Z = origin) in
#    J. More' and G. Toraldo,
#    "On the Solution of Large Quadratic-Programming Problems with Bound
#    Constraints", 
#    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
# 
#    SIF input: Ph. Toint, Dec 1989.
#    modified by Peihuang Chen, according to MINPACK-2, Apr 1992.
# 
#    classification = "QBR2-MY-V-0"
# 
#    Q is half the number of discretized points along the X axis
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TORSIOND'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'TORSIOND'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['Q'] = int(2);  #  SIF file default value
        else:
            v_['Q'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE Q                   5              $-PARAMETER n= 100
# IE Q                   11             $-PARAMETER n= 484
# IE Q                   16             $-PARAMETER n= 1024
# IE Q                   37             $-PARAMETER n= 5476
# IE Q                   50             $-PARAMETER n= 10000
# IE Q                   61             $-PARAMETER n= 14884
        if nargin<2:
            v_['C'] = float(10.0);  #  SIF file default value
        else:
            v_['C'] = float(args[1])
        v_['Q+1'] = 1+v_['Q']
        v_['P'] = v_['Q']+v_['Q']
        v_['P-1'] = -1+v_['P']
        v_['1/H'] = float(v_['P-1'])
        v_['H'] = 1.0/v_['1/H']
        v_['H2'] = v_['H']*v_['H']
        v_['C0'] = v_['H2']*v_['C']
        v_['LC'] = -1.0*v_['C0']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['P'])+1):
            for J in range(int(v_['2']),int(v_['P'])+1):
                [ig,ig_,_] = s2x_ii('GL'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2x_ii('GR'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2x_ii('G',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['LC'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            pb.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xlower[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
        for I in range(int(v_['2']),int(v_['Q'])+1):
            for J in range(int(v_['2']),int(I)+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['UPPL'] = v_['RJ-1']*v_['H']
                v_['LOWL'] = -1.0*v_['UPPL']
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWL']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPL']
            v_['MI'] = -1*I
            v_['P-I'] = v_['P']+v_['MI']
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['UPPM'] = v_['RI-1']*v_['H']
            v_['LOWM'] = -1.0*v_['UPPM']
            v_['P-I+1'] = 1+v_['P-I']
            for J in range(int(I),int(v_['P-I+1'])+1):
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWM']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPM']
            for J in range(int(v_['P-I+1']),int(v_['P-1'])+1):
                v_['MJ'] = -1*J
                v_['P-J'] = v_['P']+v_['MJ']
                v_['RP-J'] = float(v_['P-J'])
                v_['UPPR'] = v_['RP-J']*v_['H']
                v_['LOWR'] = -1.0*v_['UPPR']
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWR']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPR']
        for I in range(int(v_['Q+1']),int(v_['P-1'])+1):
            v_['MI'] = -1*I
            v_['P-I'] = v_['P']+v_['MI']
            v_['P-I+1'] = 1+v_['P-I']
            for J in range(int(v_['2']),int(v_['P-I+1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['UPPL'] = v_['RJ-1']*v_['H']
                v_['LOWL'] = -1.0*v_['UPPL']
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWL']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPL']
            v_['RP-I'] = float(v_['P-I'])
            v_['UPPM'] = v_['RP-I']*v_['H']
            v_['LOWM'] = -1.0*v_['UPPM']
            for J in range(int(v_['P-I+1']),int(I)+1):
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWM']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPM']
            for J in range(int(I),int(v_['P-1'])+1):
                v_['MJ'] = -1*J
                v_['P-J'] = v_['P']+v_['MJ']
                v_['RP-J'] = float(v_['P-J'])
                v_['UPPR'] = v_['RP-J']*v_['H']
                v_['LOWR'] = -1.0*v_['UPPR']
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOWR']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPPR']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eISQ', iet_)
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
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['P'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['2']),int(v_['P'])+1):
                v_['J-1'] = -1+J
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'D'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            for J in range(int(v_['1']),int(v_['P-1'])+1):
                ig = ig_['GR'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.25))
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.25))
        for I in range(int(v_['2']),int(v_['P'])+1):
            for J in range(int(v_['2']),int(v_['P'])+1):
                ig = ig_['GL'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.25))
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.25))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QBR2-MY-V-0"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
