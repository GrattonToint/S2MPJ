from s2mpjlib import *
class  DTOC4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC4
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 1 control variable and 2 state variables.
# 
#    The problem is not convex.
# 
#    Sources: problem 4 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    G. Di Pillo, L. Grippo and F. Lampariello,
#    "A class of structures quasi-Newton algorithms for optimal control
#    problems",
#    in H.E. Rauch, ed., IFAC Applications of nonlinear programming to
#    optimization and control, pp. 101-107, IFAC, Pergamon Press, 1983.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "QOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has 3N-1  variables (of which 2 are fixed),
#    and 2(N-1) constraints
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n=   29,m= 18 original value
# IE N                   50             $-PARAMETER  n=  149,m= 98
# IE N                   100            $-PARAMETER  n=  299,m=198
# IE N                   500            $-PARAMETER  n= 1499,m=998
# IE N                   1000           $-PARAMETER  n= 2999,m=1998
# IE N                   1500           $-PARAMETER  n= 4499,m=2998
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DTOC4'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DTOC4'
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
# IE N                   5000           $-PARAMETER  n=14999,m=9998
        v_['N-1'] = -1+v_['N']
        v_['1'] = 1
        v_['2'] = 2
        v_['RN'] = float(v_['N'])
        v_['H'] = 1.0/v_['RN']
        v_['5H'] = 5.0*v_['H']
        v_['1/5H'] = 1.0/v_['5H']
        v_['1+5H'] = 1.0+v_['5H']
        v_['-5H'] = -1.0*v_['5H']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(T),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(T))
        for T in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(T)+','+str(int(v_['1'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(T)+','+str(int(v_['1'])))
            [iv,ix_,_] = s2mpj_ii('Y'+str(T)+','+str(int(v_['2'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(T)+','+str(int(v_['2'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/5H']))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['Y'+str(T)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(v_['1+5H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['Y'+str(T)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(v_['-5H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['1'])))
            iv = ix_['X'+str(T)]
            pbm.A[ig,iv] = float(v_['5H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['2'])))
            iv = ix_['Y'+str(int(v_['T+1']))+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['Y'+str(T)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('TT'+str(T)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T)+','+str(int(v_['2'])))
            iv = ix_['Y'+str(T)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(v_['5H'])+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['2']))]] = 1.0
        pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['2']))]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = float(0.0)
        pb.x0[ix_['Y'+str(int(v_['1']))+','+str(int(v_['2']))]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2mpj_ii( 'eAAB', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'Y1SQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'Y1SQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['1']))+','+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Y2SQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'Y2SQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XSQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'XSQ'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for T in range(int(v_['2']),int(v_['N-1'])+1):
            ename = 'Y1SQ'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'Y'+str(T)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y2SQ'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'Y'+str(T)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'XSQ'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(T)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Y1SQ'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'Y1SQ'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['N']))+','+str(int(v_['1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Y2SQ'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'Y2SQ'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['N']))+','+str(int(v_['2']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'E'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eAAB')
            ielftype = arrset(ielftype, ie, iet_["eAAB"])
            vname = 'Y'+str(T)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='A')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(T)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='B')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y1SQ'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y2SQ'+str(int(v_['1']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        for T in range(int(v_['2']),int(v_['N-1'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y1SQ'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y2SQ'+str(T)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y1SQ'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y2SQ'+str(int(v_['N']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.5))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['TT'+str(T)+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-5H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  10)      3.75078392210
# LO SOLUTION(  50)      3.02963141755
# LO SOLUTION( 100)      2.94726711402
# LO SOLUTION( 500)      2.87827434035
# LO SOLUTION(1000)      2.87483889886
# LO SOLUTION(5000)      2.86386891514
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
        pb.pbclass = "QOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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
    def eAAB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]*EV_[1]
            g_[1] = EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EV_[1]+EV_[1]
                H_[0,1] = EV_[0]+EV_[0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

