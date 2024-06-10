from s2mpjlib import *
class  DTOC6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC6
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 1 control variable and 1 state variable.
# 
#    The problem is convex.
# 
#    Sources: problem 6 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    D.M. Murray and S.J. Yakowitz,
#    "The application of optimal contraol methodology to nonlinear programming
#    problems",
#    Mathematical Programming 21, pp. 331-347, 1981.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "OOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has 2N-1  variables (of which 1 is fixed),
#    and N-1 constraints
# 
#           Alternative values for the SIF file parameters:
# IE N                   11             $-PARAMETER n =   21, m =  10
# IE N                   21             $-PARAMETER n =   41, m =  20
# IE N                   31             $-PARAMETER n =   61, m =  30
# IE N                   41             $-PARAMETER n =   81, m =  40
# IE N                   51             $-PARAMETER n =  101, m =  50
# IE N                   61             $-PARAMETER n =  121, m =  60
# IE N                   71             $-PARAMETER n =  141, m =  70
# IE N                   81             $-PARAMETER n =  161, m =  80
# IE N                   91             $-PARAMETER n =  181, m =  90
# IE N                   101            $-PARAMETER n =  201, m = 100
# IE N                   501            $-PARAMETER n = 1001, m = 500
# IE N                   1001           $-PARAMETER n = 2001, m =1000
# IE N                   5001           $-PARAMETER n =10001, m =5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DTOC6'

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
            v_['N'] = int(11);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['N-1'] = -1+v_['N']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(T),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(T))
        for T in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(T),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(T))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('OY'+str(T),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['Y'+str(T)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
            [ig,ig_,_] = s2mpj_ii('OX'+str(T),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(T)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            v_['T+1'] = 1+T
            [ig,ig_,_] = s2mpj_ii('TT'+str(T),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'TT'+str(T))
            iv = ix_['Y'+str(int(v_['T+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['Y'+str(T)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['Y'+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['1']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['Y'+str(int(v_['1']))]] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'E'+str(T)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset(ielftype, ie, iet_["eEXP"])
            vname = 'X'+str(T)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for T in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['OX'+str(T)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['OY'+str(T)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['TT'+str(T)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(T)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  11)      19.80411526774
# LO SOLUTION(  21)      62.49481326823
# LO SOLUTION(  31)      119.0328455446
# LO SOLUTION(  41)      185.8961987565
# LO SOLUTION(  51)      261.1131573312
# LO SOLUTION(  61)      343.3995402405
# LO SOLUTION(  71)      431.8430677467
# LO SOLUTION(  81)      525.7575776948
# LO SOLUTION(  91)      624.6051839906
# LO SOLUTION( 101)      727.9505731659
# LO SOLUTION( 501)      6846.330143698
# LO SOLUTION(1001)      17176.03828316
# LO SOLUTION(5001)      
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
        pb.pbclass = "OOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EZ = np.exp(EV_[0])
        f_   = EZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EZ
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EZ
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

