from s2xlib import *
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
#    classification = "OOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has (N-1)*NX+N*NY  variables (of which NY are fixed),
#    and (N-1)*NY constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DTOC2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DTOC2'
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
        for T in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('O'+str(T),ig_)
            gtype = arrset(gtype,ig,'<>')
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2x_ii('TT'+str(T)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'TT'+str(T)+','+str(J))
                iv = ix_['Y'+str(int(v_['T+1']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            v_['RI'] = float(I)
            v_['TMP'] = v_['RI']*v_['1/2NY']
            pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(I)]] = v_['TMP']
            pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(I)]] = v_['TMP']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            v_['RI'] = float(I)
            v_['TMP'] = v_['RI']*v_['1/2NY']
            pb.x0[ix_['Y'+str(int(v_['1']))+','+str(I)]] = float(v_['TMP'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eOEL', iet_)
        elftv = loaset(elftv,it,0,'YY1')
        elftv = loaset(elftv,it,1,'YY2')
        elftv = loaset(elftv,it,2,'YY3')
        elftv = loaset(elftv,it,3,'YY4')
        elftv = loaset(elftv,it,4,'XX1')
        elftv = loaset(elftv,it,5,'XX2')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'YY')
        [it,iet_,_] = s2x_ii( 'eSINE', iet_)
        elftv = loaset(elftv,it,0,'ZZ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'EO'+str(T)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eOEL')
            ielftype = arrset(ielftype, ie, iet_["eOEL"])
            vname = 'Y'+str(T)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='YY1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='YY2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['3']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='YY3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['4']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='YY4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(T)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XX1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(T)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XX2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'SY'+str(T)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSINE')
                ielftype = arrset(ielftype, ie, iet_["eSINE"])
                vname = 'Y'+str(T)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='ZZ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ename = 'SX'+str(T)+','+str(I)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSINE')
                ielftype = arrset(ielftype, ie, iet_["eSINE"])
                vname = 'X'+str(T)+','+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='ZZ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['NY'])+1):
            ename = 'YNSQ'+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'Y'+str(int(v_['N']))+','+str(J)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['O'+str(T)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EO'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['NY'])+1):
            ig = ig_['O'+str(int(v_['N']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['YNSQ'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['TT'+str(T)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SY'+str(T)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                for I in range(int(v_['1']),int(v_['NX'])+1):
                    ig = ig_['TT'+str(T)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SX'+str(T)+','+str(I)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['C'+str(J)+','+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
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
        pb.pbclass = "OOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
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
    def eSINE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SZ = np.sin(EV_[0])
        f_   = SZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(EV_[0])
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
    def eOEL(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XN2 = EV_[4]*EV_[4]+EV_[5]*EV_[5]
        YN2 = EV_[0]*EV_[0]+EV_[1]*EV_[1]+EV_[2]*EV_[2]+EV_[3]*EV_[3]
        SZ = np.sin(0.5*XN2)
        CZ = np.cos(0.5*XN2)
        SZ2 = SZ*SZ+1.0
        SC = SZ*CZ
        CCSS = CZ*CZ-SZ*SZ
        f_   = YN2*SZ2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[4] = 2.0*YN2*SC*EV_[4]
            g_[5] = 2.0*YN2*SC*EV_[5]
            g_[0] = 2.0*EV_[0]*SZ2
            g_[1] = 2.0*EV_[1]*SZ2
            g_[2] = 2.0*EV_[2]*SZ2
            g_[3] = 2.0*EV_[3]*SZ2
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[4,4] = 2.0*YN2*(SC+EV_[4]*EV_[4]*CCSS)
                H_[4,5] = 2.0*YN2*EV_[4]*EV_[5]*CCSS
                H_[5,4] = H_[4,5]
                H_[4,0] = 4.0*EV_[0]*SC*EV_[4]
                H_[0,4] = H_[4,0]
                H_[4,1] = 4.0*EV_[1]*SC*EV_[4]
                H_[1,4] = H_[4,1]
                H_[4,2] = 4.0*EV_[2]*SC*EV_[4]
                H_[2,4] = H_[4,2]
                H_[4,3] = 4.0*EV_[3]*SC*EV_[4]
                H_[3,4] = H_[4,3]
                H_[5,5] = 2.0*YN2*(SC+EV_[5]*EV_[5]*CCSS)
                H_[5,0] = 4.0*EV_[0]*SC*EV_[5]
                H_[0,5] = H_[5,0]
                H_[5,1] = 4.0*EV_[1]*SC*EV_[5]
                H_[1,5] = H_[5,1]
                H_[5,2] = 4.0*EV_[2]*SC*EV_[5]
                H_[2,5] = H_[5,2]
                H_[5,3] = 4.0*EV_[3]*SC*EV_[5]
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

