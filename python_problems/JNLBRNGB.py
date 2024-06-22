from s2mpjlib import *
class  JNLBRNGB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : JNLBRNGB
#    *********
# 
#    The quadratic journal bearing problem (with excentricity = 0.5)
# 
#    Source:
#    J. More' and G. Toraldo,
#    "On the Solution of Large Quadratic-Programming Problems with Bound
#    Constraints", 
#    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "QBR2-AY-V-0"
# 
#    The rectangle is discretized into (pt-1)(py-1) little rectangles. The
#    heights of the considered surface above the corners of these little
#    rectangles are the problem variables,  There are px*py of them.
# 
#    PT is the number of points along the T (\theta) side of the rectangle
#    PY is the number of points along the Y side of the rectangle
# 
#           Alternative values for the SIF file parameters:
# IE PT                  4              $-PARAMETER  n=16
# IE PY                  4              $-PARAMETER
# 
# IE PT                  10             $-PARAMETER  n=100
# IE PY                  10             $-PARAMETER
# 
# IE PT                  23             $-PARAMETER  n=529
# IE PY                  23             $-PARAMETER
# 
# IE PT                  32             $-PARAMETER  n=1024
# IE PY                  32             $-PARAMETER
# 
# IE PT                  34             $-PARAMETER  n=1156
# IE PY                  34             $-PARAMETER
# 
# IE PT                  75             $-PARAMETER  n=5625   original value
# IE PY                  75             $-PARAMETER           original value
# 
# IE PT                  100            $-PARAMETER  n=10000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'JNLBRNGB'

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
            v_['PT'] = int(5);  #  SIF file default value
        else:
            v_['PT'] = int(args[0])
# IE PY                  100            $-PARAMETER
        if nargin<2:
            v_['PY'] = int(5);  #  SIF file default value
        else:
            v_['PY'] = int(args[1])
# IE PT                  125            $-PARAMETER  n=15625
# IE PY                  125            $-PARAMETER
        if nargin<3:
            v_['EX'] = float(0.5);  #  SIF file default value
        else:
            v_['EX'] = float(args[2])
        v_['LT'] = 6.2831853
        v_['LY'] = 20.0
        v_['PT-1'] = -1+v_['PT']
        v_['RPT-1'] = float(v_['PT-1'])
        v_['HT1'] = 1.0/v_['RPT-1']
        v_['HT'] = v_['HT1']*v_['LT']
        v_['1/HT'] = 1.0/v_['HT']
        v_['PY-1'] = -1+v_['PY']
        v_['RPY-1'] = float(v_['PY-1'])
        v_['HY1'] = 1.0/v_['RPY-1']
        v_['HY'] = v_['HY1']*v_['LY']
        v_['1/HY'] = 1.0/v_['HY']
        v_['HTHY'] = v_['HT']*v_['HY']
        v_['HT/HY'] = v_['HT']*v_['1/HY']
        v_['HY/HT'] = v_['HY']*v_['1/HT']
        v_['EXHTHY'] = v_['HTHY']*v_['EX']
        v_['CLINC'] = -1.0*v_['EXHTHY']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['PT'])+1):
            for J in range(int(v_['1']),int(v_['PY'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['PT-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XI1'] = v_['RI-1']*v_['HT']
            v_['SXI1'] = np.sin(v_['XI1'])
            v_['COEFF'] = v_['SXI1']*v_['CLINC']
            for J in range(int(v_['2']),int(v_['PY-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['COEFF'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['PY'])+1):
            pb.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xlower[ix_['X'+str(int(v_['PT']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['PT']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['PT-1'])+1):
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['PY']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['PY']))]] = 0.0
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
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
        for I in range(int(v_['2']),int(v_['PT-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            for J in range(int(v_['2']),int(v_['PY-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'D'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['PT-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XI1'] = v_['RI-1']*v_['HT']
            v_['CXI1'] = np.cos(v_['XI1'])
            v_['ECX'] = v_['CXI1']*v_['EX']
            v_['ECX1'] = 1.0+v_['ECX']
            v_['E12'] = v_['ECX1']*v_['ECX1']
            v_['WI'] = v_['ECX1']*v_['E12']
            v_['2WI'] = v_['WI']+v_['WI']
            v_['RI'] = float(I)
            v_['XI+1'] = v_['RI']*v_['HT']
            v_['CXI+1'] = np.cos(v_['XI+1'])
            v_['E+CX0'] = v_['CXI+1']*v_['EX']
            v_['E+CX1'] = 1.0+v_['E+CX0']
            v_['E22'] = v_['E+CX1']*v_['E+CX1']
            v_['WI+1'] = v_['E+CX1']*v_['E22']
            v_['I-2'] = -2+I
            v_['RI-2'] = float(v_['I-2'])
            v_['XI-1'] = v_['RI-2']*v_['HT']
            v_['CXI-1'] = np.cos(v_['XI-1'])
            v_['E-CX0'] = v_['CXI-1']*v_['EX']
            v_['E-CX1'] = 1.0+v_['E-CX0']
            v_['E32'] = v_['E-CX1']*v_['E-CX1']
            v_['WI-1'] = v_['E-CX1']*v_['E32']
            v_['PM0'] = v_['2WI']*v_['WI+1']
            v_['PM1'] = 0.0833333333*v_['PM0']
            v_['MU/HY2'] = v_['PM1']*v_['HT/HY']
            v_['MU/HT2'] = v_['PM1']*v_['HY/HT']
            v_['PL0'] = v_['2WI']*v_['WI-1']
            v_['PL1'] = 0.0833333333*v_['PL0']
            v_['LA/HY2'] = v_['PL1']*v_['HT/HY']
            v_['LA/HT2'] = v_['PL1']*v_['HY/HT']
            for J in range(int(v_['2']),int(v_['PY-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MU/HT2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['MU/HY2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['LA/HT2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['LA/HY2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(4)            -1.8551D+01
# LO SOLTN(10)           -7.2552D+00
# LO SOLTN(23)           -6.5068D+00
# LO SOLTN(32)           -6.4401D+00
# LO SOLTN(75)           -6.3297D+00
# LO SOLTN(100)          -6.3007D+00
# LO SOLTN(125)          -6.2807D+00
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QBR2-AY-V-0"
        pb.x0          = np.zeros((pb.n,1))
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

