from s2xlib import *
class  DTOC1NA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC1NA
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, NX control variables and NY state variables.
#    The nonlinearity parameter mu is set to 0.005.
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
#    classification = "OQR2-AN-V-V"
# 
#    Problem variants: they are identified by the values of
#    the parameter vector ( N, NX, NY )
# 
#    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
#    and (N-1)*NY constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DTOC1NA'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DTOC1NA'
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
            v_['NX'] = int(2);  #  SIF file default value
        else:
            v_['NX'] = int(args[1])
        if nargin<3:
            v_['NY'] = int(4);  #  SIF file default value
        else:
            v_['NY'] = int(args[2])
#           Alternative values for the SIF file parameters:
# IE N                   50             $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  298, m= 196
# IE NY                  4              $-PARAMETER # states   }
# IE N                   100            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n=  598, m= 396
# IE NY                  4              $-PARAMETER # states   }
# IE N                   500            $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n= 2998, m=1996
# IE NY                  4              $-PARAMETER # states   }
# IE N                   1000           $-PARAMETER # periods  }
# IE NX                  2              $-PARAMETER # controls } n= 5998, m=3996
# IE NY                  4              $-PARAMETER # states   }
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
            v_['MU'] = float(0.005);  #  SIF file default value
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(T)+','+str(I),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(T)+','+str(I))
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2x_ii('Y'+str(T)+','+str(I),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(T)+','+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2x_ii('OX'+str(T)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(T)+','+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2x_ii('OY'+str(T)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Y'+str(T)+','+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['Y'+str(T)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
            iv = ix_['Y'+str(T)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(0.25)+pbm.A[ig,iv]
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
                iv = ix_['X'+str(T)+','+str(I)]
                pbm.A[ig,iv] = float(v_['B'+str(int(v_['1']))+','+str(I)])+pbm.A[ig,iv]
            for J in range(int(v_['2']),int(v_['NY-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                iv = ix_['Y'+str(int(v_['T+1']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(T)+','+str(J)]
                pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
                iv = ix_['Y'+str(T)+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(-0.25)+pbm.A[ig,iv]
                iv = ix_['Y'+str(T)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(0.25)+pbm.A[ig,iv]
                for I in range(int(v_['1']),int(v_['NX'])+1):
                    [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                    iv = ix_['X'+str(T)+','+str(I)]
                    pbm.A[ig,iv] = float(v_['B'+str(J)+','+str(I)])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
            iv = ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['NY']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
            iv = ix_['Y'+str(T)+','+str(int(v_['NY']))]
            pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
            iv = ix_['Y'+str(T)+','+str(int(v_['NY-1']))]
            pbm.A[ig,iv] = float(-0.25)+pbm.A[ig,iv]
            for I in range(int(v_['1']),int(v_['NX'])+1):
                [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(int(v_['NY'])),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['NY'])))
                iv = ix_['X'+str(T)+','+str(I)]
                pbm.A[ig,iv] = float(v_['B'+str(int(v_['NY']))+','+str(I)])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['OX'+str(T)+','+str(I)],float(-0.5))
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['OY'+str(T)+','+str(I)],float(-0.25))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(I)]] = 0.0
            pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'ePR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'MUC')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for K in range(int(v_['0']),int(v_['NYNX-1'])+1):
                v_['I'] = int(np.fix(K/v_['NX']))
                v_['INX'] = v_['I']*v_['NX']
                v_['J'] = K-v_['INX']
                v_['I'] = 1+v_['I']
                v_['J'] = 1+v_['J']
                ename = 'E'+str(T)+','+str(K)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePR')
                ielftype = arrset(ielftype, ie, iet_["ePR"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'Y'+str(T)+','+str(int(v_['I']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(T)+','+str(int(v_['J']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='MUC')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(int(v_['I']))+','+str(int(v_['J']))])))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL4',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ig = ig_['OX'+str(T)+','+str(I)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL4')
        for T in range(int(v_['1']),int(v_['N'])+1):
            for I in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['OY'+str(T)+','+str(I)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL4')
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                for K in range(int(v_['0']),int(v_['NYNX-1'])+1):
                    ig = ig_['TT'+str(T)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(T)+','+str(K)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OQR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EV_[1]
            g_[1] = pbm.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL4(pbm,nargout,*args):

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

