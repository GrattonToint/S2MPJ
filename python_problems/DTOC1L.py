from s2mpjlib import *
class  DTOC1L(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC1L
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, NX control variables and NY state variables.
#    The parameter mu in the original problem formulation is set to zero, 
#    yielding linear transition functions, hence the L in the problem's name.
# 
#    The problem is convex.
# 
#    Sources: problem 1 (with mu = 0) in
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
#    classification = "C-COLR2-AN-V-V"
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

    name = 'DTOC1L'

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
        v_['NY-1'] = -1+v_['NY']
        v_['NX+NY'] = v_['NX']+v_['NY']
        v_['RXY'] = float(v_['NX+NY'])
        v_['1/RXY'] = 1.0/v_['RXY']
        for J in range(int(v_['1']),int(v_['NX'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                v_['I-J'] = I-J
                v_['RI-J'] = float(v_['I-J'])
                v_['B'+str(I)+','+str(J)] = v_['RI-J']*v_['1/RXY']
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
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
# LO SOLUTION(  10,2, 4) 0.0735931360
# LO SOLUTION(  50,2, 4) 0.2299411960
# LO SOLUTION( 100,2, 4) 0.4253681120
# LO SOLUTION( 500,2, 4) 1.9887794988
# LO SOLUTION(1000,2, 4) 3.9430507151
# LO SOLUTION(  10,5,10) 1.1498579294
# LO SOLUTION(  50,5,10) 6.1678479713
# LO SOLUTION( 100,5,10) 12.439954329
# LO SOLUTION( 500,5,10) 62.616843379
# LO SOLUTION(1000,5,10) 125.33793359
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-COLR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

