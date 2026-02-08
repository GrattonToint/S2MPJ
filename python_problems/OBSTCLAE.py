from s2mpjlib import *
class  OBSTCLAE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OBSTCLAE
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
#    See also More 1989 (Problem A, Starting point E)
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CQBR2-AY-V-0"
# 
#    PX is the number of points along the X side of the rectangle
#    PY is the number of points along the Y side of the rectangle
# 
#           Alternative values for the SIF file parameters:
# IE PX                  4              $-PARAMETER n = 16
# IE PY                  4              $-PARAMETER
# 
# IE PX                  10             $-PARAMETER n = 100
# IE PY                  10             $-PARAMETER
# 
# IE PX                  23             $-PARAMETER n = 529
# IE PY                  23             $-PARAMETER
# 
# IE PX                  32             $-PARAMETER n = 1024
# IE PY                  32             $-PARAMETER
# 
# IE PX                  75             $-PARAMETER n = 5625     original value
# IE PY                  75             $-PARAMETER              original value
# 
# IE PX                  100            $-PARAMETER n = 10000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OBSTCLAE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['PX'] = int(5);  #  SIF file default value
        else:
            v_['PX'] = int(args[0])
# IE PY                  100            $-PARAMETER
        if nargin<2:
            v_['PY'] = int(20);  #  SIF file default value
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
        v_['LINC'] = -1.0*v_['C0']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['PX'])+1):
            for I in range(int(v_['1']),int(v_['PY'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['LINC']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),2000.0)
        self.xlower = np.zeros((self.n,1))
        for J in range(int(v_['1']),int(v_['PX'])+1):
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xlower[ix_['X'+str(int(v_['PY']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['PY']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            self.xlower[ix_['X'+str(I)+','+str(int(v_['PX']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['PX']))]] = 0.0
            self.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['XI1'] = v_['RI-1']*v_['HY']
            v_['3XI1'] = 3.2*v_['XI1']
            v_['SXI1'] = np.sin(v_['3XI1'])
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['XI2'] = v_['RJ-1']*v_['HX']
                v_['3XI2'] = 3.3*v_['XI2']
                v_['SXI2'] = np.sin(v_['3XI2'])
                v_['LOW'] = v_['SXI1']*v_['SXI2']
                self.xlower[ix_['X'+str(I)+','+str(J)]] = v_['LOW']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        for J in range(int(v_['1']),int(v_['PX'])+1):
            self.x0[ix_['X'+str(int(v_['1']))+','+str(J)]] = float(0.0)
            self.x0[ix_['X'+str(int(v_['PY']))+','+str(J)]] = float(0.0)
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            self.x0[ix_['X'+str(I)+','+str(int(v_['PX']))]] = float(0.0)
            self.x0[ix_['X'+str(I)+','+str(int(v_['1']))]] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'D'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eISQ')
                ielftype = arrset(ielftype,ie,iet_["eISQ"])
                vname = 'X'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(2000.0),float(1.0))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['PY-1'])+1):
            for J in range(int(v_['2']),int(v_['PX-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['HY/4HX']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['HX/4HY']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['HY/4HX']))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel,float(v_['HX/4HY']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)            0.753659753
# LO SOLTN(10)           1.397897560
# LO SOLTN(23)           1.678027027
# LO SOLTN(32)           1.748270031
# LO SOLTN(75)           ???
# LO SOLTN(100)          ???
# LO SOLTN(125)          ???
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CQBR2-AY-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        f_   = IV_[0]*IV_[0]
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

