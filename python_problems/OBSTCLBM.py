from s2mpjlib import *
class  OBSTCLBM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OBSTCLBM
#    *********
# 
#    A quadratic obstacle problem by Dembo and Tulowitzki
# 
#    The problem comes from the obstacle problem on a rectangle.
#    The rectangle is discretized into (px-1)(py-1) little rectangles. The
#    heights of the considered surface above the corners of these little
#    rectangles are the problem variables,  There are px*py of them.
# 
#    Source:
#    R. Dembo and U. Tulowitzki,
#    "On the minimization of quadratic functions subject to box
#    constraints",
#    WP 71, Yale University (new Haven, USA), 1983.
# 
#    See also More 1989 (Problem B, Starting point M (average of L and U))
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "QBR2-AY-V-0"
# 
#    PX is the number of points along the X side of the rectangle
#    PY is the number of points along the Y side of the rectangle
# 
#           Alternative values for the SIF file parameters:
# IE PX                  4              $-PARAMETER n = 16
# IE PY                  4              $-PARAMETER
# 
# IE PX                  10             $-PARAMETER n = 100     original value
# IE PY                  10             $-PARAMETER             original value
# 
# IE PX                  23             $-PARAMETER n = 529
# IE PY                  23             $-PARAMETER
# 
# IE PX                  32             $-PARAMETER n = 1024
# IE PY                  32             $-PARAMETER
# 
# IE PX                  75             $-PARAMETER n = 5625
# IE PY                  75             $-PARAMETER
# 
# IE PX                  100            $-PARAMETER n = 10000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OBSTCLBM'

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
            v_['PX'] = int(5);  #  SIF file default value
        else:
            v_['PX'] = int(args[0])
        if nargin<2:
            v_['PY'] = int(100);  #  SIF file default value
        else:
            v_['PY'] = int(args[1])
# IE PX                  125            $-PARAMETER n = 15625
# IE PY                  125            $-PARAMETER
        if nargin<3:
            v_['C'] = float(1.0);  #  SIF file default value
        else:
            v_['C'] = float(args[2])
        v_['PX-1'] = -1+v_['PX']
        v_['RPX-1'] = float(v_['PX-1'])
        v_['HX'] = 1.0/v_['RPX-1']
        v_['PY-1'] = -1+v_['PY']
        v_['RPY-1'] = float(v_['PY-1'])
        v_['HY'] = 1.0/v_['RPY-1']
        v_['HXHY'] = v_['HX']*v_['HY']
        v_['1/HX'] = 1.0/v_['HX']
        v_['1/HY'] = 1.0/v_['HY']
        v_['HX/HY'] = v_['HX']*v_['1/HY']
        v_['HY/HX'] = v_['HY']*v_['1/HX']
        v_['HY/4HX'] = 0.25*v_['HY/HX']
        v_['HX/4HY'] = 0.25*v_['HX/HY']
        v_['C0'] = v_['HXHY']*v_['C']
        v_['LC'] = -1.0*v_['C0']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['PX'])+1):
            for I in range(int(v_['1']),int(v_['PY'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['LC'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['PX'])+1):
            pb.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xlower[ix_['X'+str(int(v_['PY']))+','+str(J)]] = 0.0
            pb.xupper[ix_['X'+str(int(v_['PY']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['PX']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['PX']))]] = 0.0
            pb.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XSI1'] = v_['RI-1']*v_['HY']
            v_['3XSI1'] = 9.2*v_['XSI1']
            v_['SXSI1'] = np.sin(v_['3XSI1'])
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['XSI2'] = v_['RJ-1']*v_['HX']
                v_['3XSI2'] = 9.3*v_['XSI2']
                v_['SXSI2'] = np.sin(v_['3XSI2'])
                v_['L1'] = v_['SXSI1']*v_['SXSI2']
                v_['L2'] = v_['L1']*v_['L1']
                v_['LOW'] = v_['L2']*v_['L1']
                pb.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOW']
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XSI1'] = v_['RI-1']*v_['HY']
            v_['3XSI1'] = 9.2*v_['XSI1']
            v_['SXSI1'] = np.sin(v_['3XSI1'])
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['XSI2'] = v_['RJ-1']*v_['HX']
                v_['3XSI2'] = 9.3*v_['XSI2']
                v_['SXSI2'] = np.sin(v_['3XSI2'])
                v_['L1'] = v_['SXSI1']*v_['SXSI2']
                v_['L2'] = v_['L1']*v_['L1']
                v_['UPP'] = 0.02+v_['L2']
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = v_['UPP']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for J in range(int(v_['1']),int(v_['PX'])+1):
            pb.x0[ix_['X'+str(int(v_['1']))+','+str(J)]] = float(0.0)
            pb.x0[ix_['X'+str(int(v_['PY']))+','+str(J)]] = float(0.0)
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            pb.x0[ix_['X'+str(I)+','+str(int(v_['PX']))]] = float(0.0)
            pb.x0[ix_['X'+str(I)+','+str(int(v_['1']))]] = float(0.0)
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XSI1'] = v_['RI-1']*v_['HY']
            v_['3XSI1'] = 9.2*v_['XSI1']
            v_['SXSI1'] = np.sin(v_['3XSI1'])
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['XSI2'] = v_['RJ-1']*v_['HX']
                v_['3XSI2'] = 9.3*v_['XSI2']
                v_['SXSI2'] = np.sin(v_['3XSI2'])
                v_['L1'] = v_['SXSI1']*v_['SXSI2']
                v_['L2'] = v_['L1']*v_['L1']
                v_['LOW'] = v_['L1']*v_['L2']
                v_['UPP'] = 0.02+v_['L2']
                v_['M1'] = v_['LOW']+v_['UPP']
                v_['MID'] = 0.5*v_['M1']
                pb.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['MID'])
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
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
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
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['HY/4HX']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['HX/4HY']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['HY/4HX']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['HX/4HY']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)            -0.0081108
# LO SOLTN(10)           2.87503823
# LO SOLTN(23)           6.51932527
# LO SOLTN(32)           6.88708670
# LO SOLTN(75)           ???
# LO SOLTN(100)          ???
# LO SOLTN(125)          ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QBR2-AY-V-0"
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

