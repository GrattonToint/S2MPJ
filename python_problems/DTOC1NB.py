from s2mpjlib import *
class  DTOC1NB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC1NB
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, NX control variables and NY state variables.
#    The nonlinearity parameter mu is set to 0.05.
# 
#    The problem is not convex.
# 
#    Sources: problem 1 in
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
#    classification = "C-COQR2-AN-V-V"
# 
#    Problem variants: they are identified by the values of
#    the parameter vector ( N, NX, NY )
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

    name = 'DTOC1NB'

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
        if nargin<4:
            v_['MU'] = float(0.05);  #  SIF file default value
        else:
            v_['MU'] = float(args[3])
        v_['N-1'] = -1+v_['N']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['NY-1'] = -1+v_['NY']
        v_['NX+NY'] = v_['NX']+v_['NY']
        v_['RXY'] = float(v_['NX+NY'])
        v_['1/RXY'] = 1.0/v_['RXY']
        v_['MU/RXY'] = v_['MU']*v_['1/RXY']
        v_['NYNX'] = v_['NX']*v_['NY']
        v_['NYNX-1'] = -1+v_['NYNX']
        for J in range(int(v_['1']),int(v_['NX'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                v_['I-J'] = I-J
                v_['RI-J'] = float(v_['I-J'])
                v_['B'+str(I)+','+str(J)] = v_['RI-J']*v_['1/RXY']
                v_['I+J'] = I+J
                v_['RI+J'] = float(v_['I+J'])
                v_['C'+str(I)+','+str(J)] = v_['RI+J']*v_['MU/RXY']
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
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2mpj_ii('OX'+str(T)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(T)+','+str(I)]])
                valA = np.append(valA,float(1.0))
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('OY'+str(T)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(T)+','+str(I)]])
                valA = np.append(valA,float(1.0))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['1']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['1']))]])
            valA = np.append(valA,float(0.5))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(0.25))
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(T)+','+str(I)]])
                valA = np.append(valA,float(v_['B'+str(int(v_['1']))+','+str(I)]))
            for J in range(int(v_['2']),int(v_['NY-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(int(v_['T+1']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(T)+','+str(J)]])
                valA = np.append(valA,float(0.5))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['J-1']))]])
                valA = np.append(valA,float(-0.25))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['J+1']))]])
                valA = np.append(valA,float(0.25))
                for I in range(int(v_['1']),int(v_['NX'])+1):
                    [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(T)+','+str(I)]])
                    valA = np.append(valA,float(v_['B'+str(J)+','+str(I)]))
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['NY']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['NY']))]])
            valA = np.append(valA,float(0.5))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(T)+','+str(int(v_['NY-1']))]])
            valA = np.append(valA,float(-0.25))
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(T)+','+str(I)]])
                valA = np.append(valA,float(v_['B'+str(int(v_['NY']))+','+str(I)]))
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                self.gconst = arrset(self.gconst,ig_['OX'+str(T)+','+str(I)],float(-0.5))
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                self.gconst = arrset(self.gconst,ig_['OY'+str(T)+','+str(I)],float(-0.25))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            self.xlower[ix_['Y'+str(int(v_['1']))+','+str(I)]] = 0.0
            self.xupper[ix_['Y'+str(int(v_['1']))+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'MUC')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for K in range(int(v_['0']),int(v_['NYNX-1'])+1):
                v_['I'] = int(np.fix(K/v_['NX']))
                v_['INX'] = v_['I']*v_['NX']
                v_['J'] = K-v_['INX']
                v_['I'] = 1+v_['I']
                v_['J'] = 1+v_['J']
                ename = 'E'+str(T)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePR')
                ielftype = arrset(ielftype,ie,iet_["ePR"])
                self.x0 = np.zeros((self.n,1))
                vname = 'Y'+str(T)+','+str(int(v_['I']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(T)+','+str(int(v_['J']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='MUC')[0]
                self.elpar  = (
                      loaset(self.elpar,ie,posep[0],float(v_['C'+str(int(v_['I']))+','+str(int(v_['J']))])))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL4',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ig = ig_['OX'+str(T)+','+str(I)]
                self.grftype = arrset(self.grftype,ig,'gL4')
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['OY'+str(T)+','+str(I)]
                self.grftype = arrset(self.grftype,ig,'gL4')
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                for K in range(int(v_['0']),int(v_['NYNX-1'])+1):
                    ig = ig_['TT'+str(T)+','+str(J)]
                    posel = len(self.grelt[ig])
                    self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(T)+','+str(K)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
# LO S(  10,2, 4)        0.09848124993
# LO S(  50,2, 4)        0.38293273441
# LO S( 100,2, 4)        0.73850718221
# LO S( 500,2, 4)        3.58309876144
# LO S(1000,2, 4)        7.13884991815
# LO S(  10,5,10)        1.40183631317
# LO S(  50,5,10)        7.86232111363
# LO S( 100,5,10)        15.9377918009
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
        self.pbclass   = "C-COQR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*EV_[1,0]
            g_[1] = self.elpar[iel_][0]*EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = self.elpar[iel_][0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL4(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**4
        if nargout>1:
            g_ = 4.0*GVAR_**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 12.0*GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

