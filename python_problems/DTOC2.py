from s2mpjlib import *
class  DTOC2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC2
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 2 control variables and 4 state variables.
# 
#    The problem is not convex.
# 
#    Sources: problem 2 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    L.Z. Liao and C.A. Shoemaker,
#    "Advantages of differential dynamic programming over Newton's method for
#    discrete-time optimal control problems",
#    Tech. Report ctc92tr97, Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "C-COOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
#    and (N-1)*NY constraints
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER # periods  } original value
# IE NX                  2              $-PARAMETER # controls } n=   58, m=  36
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   50             $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  298, m= 196
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   100            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  598, m= 396
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   500            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n= 2998, m=1996
# IE NY                  4              $-PARAMETER # states   }
# 
# IE N                   1000           $-PARAMETER # periods  }
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DTOC2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE NX                  2              $-PARAMETER # controls } n= 5998, m=3996
        if nargin<2:
            v_['NX'] = int(2);  #  SIF file default value
        else:
            v_['NX'] = int(args[1])
# IE NY                  4              $-PARAMETER # states   }
        if nargin<3:
            v_['NY'] = int(4);  #  SIF file default value
        else:
            v_['NY'] = int(args[2])
# IE N                   10             $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=  145, m=  90
# IE NY                  10             $-PARAMETER # states   }
# IE N                   50             $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=  745, m= 490
# IE NY                  10             $-PARAMETER # states   }
# IE N                   100            $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n= 1495, m= 990
# IE NY                  10             $-PARAMETER # states   }
# IE N                   500            $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n= 7495, m=4990
# IE NY                  10             $-PARAMETER # states   }
# IE N                   1000           $-PARAMETER # periods  }
# IE NX                  5              $-PARAMETER # controls } n=14995, m=9990
# IE NY                  10             $-PARAMETER # states   }
        v_['N-1'] = -1+v_['N']
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['NY-1'] = -1+v_['NY']
        v_['2NY'] = v_['NY']+v_['NY']
        v_['R2NY'] = float(v_['2NY'])
        v_['1/2NY'] = 1.0/v_['R2NY']
        for J in range(int(v_['1']),int(v_['NX'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                v_['I+J'] = I+J
                v_['RI+J'] = float(v_['I+J'])
                v_['C'+str(I)+','+str(J)] = v_['RI+J']*v_['1/2NY']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(T)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(T)+','+str(I))
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('Y'+str(T)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'Y'+str(T)+','+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for T in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(T),ig_)
            gtype = arrset(gtype,ig,'<>')
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['T+1']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            v_['RI'] = float(I)
            v_['TMP'] = v_['RI']*v_['1/2NY']
            self.xlower[ix_['Y'+str(int(v_['1']))+','+str(I)]] = v_['TMP']
            self.xupper[ix_['Y'+str(int(v_['1']))+','+str(I)]] = v_['TMP']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            v_['RI'] = float(I)
            v_['TMP'] = v_['RI']*v_['1/2NY']
            self.x0[ix_['Y'+str(int(v_['1']))+','+str(I)]] = float(v_['TMP'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOEL', iet_)
        elftv = loaset(elftv,it,0,'YY1')
        elftv = loaset(elftv,it,1,'YY2')
        elftv = loaset(elftv,it,2,'YY3')
        elftv = loaset(elftv,it,3,'YY4')
        elftv = loaset(elftv,it,4,'XX1')
        elftv = loaset(elftv,it,5,'XX2')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'YY')
        [it,iet_,_] = s2mpj_ii( 'eSINE', iet_)
        elftv = loaset(elftv,it,0,'ZZ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'EO'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eOEL')
            ielftype = arrset(ielftype,ie,iet_["eOEL"])
            vname = 'Y'+str(T)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='YY1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='YY2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='YY3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['4']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='YY4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(T)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XX1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(T)+','+str(int(v_['2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XX2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'SY'+str(T)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSINE')
                ielftype = arrset(ielftype,ie,iet_["eSINE"])
                vname = 'Y'+str(T)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='ZZ')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ename = 'SX'+str(T)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSINE')
                ielftype = arrset(ielftype,ie,iet_["eSINE"])
                vname = 'X'+str(T)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='ZZ')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['NY'])+1):
            ename = 'YNSQ'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'Y'+str(int(v_['N']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='YY')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['O'+str(T)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EO'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['NY'])+1):
            ig = ig_['O'+str(int(v_['N']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['YNSQ'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['TT'+str(T)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['SY'+str(T)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                for I in range(int(v_['1']),int(v_['NX'])+1):
                    ig = ig_['TT'+str(T)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt = loaset(self.grelt,ig,posel,ie_['SX'+str(T)+','+str(I)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,float(v_['C'+str(J)+','+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
# LO SOLUTION(  10)      0.485983918948
# LO SOLUTION(  20)      0.486212213803
# LO SOLUTION(  30)      0.486383392574
# LO SOLUTION(  40)      0.486572686778
# LO SOLUTION(  50)      0.486884900389
# LO SOLUTION( 100)      0.487532342563
# LO SOLUTION( 500)      0.490996540460
# LO SOLUTION(1000)      0.490200910983
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSINE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SZ = np.sin(EV_[0,0])
        f_   = SZ
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -SZ
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOEL(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XN2 = EV_[4,0]*EV_[4,0]+EV_[5,0]*EV_[5,0]
        YN2  = (
              EV_[0,0]*EV_[0,0]+EV_[1,0]*EV_[1,0]+EV_[2,0]*EV_[2,0]+EV_[3,0]*EV_[3,0])
        SZ = np.sin(0.5*XN2)
        CZ = np.cos(0.5*XN2)
        SZ2 = SZ*SZ+1.0
        SC = SZ*CZ
        CCSS = CZ*CZ-SZ*SZ
        f_   = YN2*SZ2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[4] = 2.0*YN2*SC*EV_[4,0]
            g_[5] = 2.0*YN2*SC*EV_[5,0]
            g_[0] = 2.0*EV_[0,0]*SZ2
            g_[1] = 2.0*EV_[1,0]*SZ2
            g_[2] = 2.0*EV_[2,0]*SZ2
            g_[3] = 2.0*EV_[3,0]*SZ2
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[4,4] = 2.0*YN2*(SC+EV_[4,0]*EV_[4,0]*CCSS)
                H_[4,5] = 2.0*YN2*EV_[4,0]*EV_[5,0]*CCSS
                H_[5,4] = H_[4,5]
                H_[4,0] = 4.0*EV_[0,0]*SC*EV_[4,0]
                H_[0,4] = H_[4,0]
                H_[4,1] = 4.0*EV_[1,0]*SC*EV_[4,0]
                H_[1,4] = H_[4,1]
                H_[4,2] = 4.0*EV_[2,0]*SC*EV_[4,0]
                H_[2,4] = H_[4,2]
                H_[4,3] = 4.0*EV_[3,0]*SC*EV_[4,0]
                H_[3,4] = H_[4,3]
                H_[5,5] = 2.0*YN2*(SC+EV_[5,0]*EV_[5,0]*CCSS)
                H_[5,0] = 4.0*EV_[0,0]*SC*EV_[5,0]
                H_[0,5] = H_[5,0]
                H_[5,1] = 4.0*EV_[1,0]*SC*EV_[5,0]
                H_[1,5] = H_[5,1]
                H_[5,2] = 4.0*EV_[2,0]*SC*EV_[5,0]
                H_[2,5] = H_[5,2]
                H_[5,3] = 4.0*EV_[3,0]*SC*EV_[5,0]
                H_[3,5] = H_[5,3]
                H_[0,0] = 2.0*SZ2
                H_[1,1] = 2.0*SZ2
                H_[2,2] = 2.0*SZ2
                H_[3,3] = 2.0*SZ2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

