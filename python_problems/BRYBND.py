from s2xlib import *
class  BRYBND(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BRYBND
#    *********
#    Broyden banded system of nonlinear equations, considered in the
#    least square sense.
# 
#    Source: problem 31 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#73 (p. 41) and Toint#18
# 
#    SDIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-AN-V-0"
# 
#    N is the number of equations and variables (variable).
# 
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BRYBND'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'BRYBND'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5000);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        if nargin<2:
            v_['KAPPA1'] = 2.0;  #  SIF file default value
        else:
            v_['KAPPA1'] = args[1]
        if nargin<3:
            v_['KAPPA2'] = 5.0;  #  SIF file default value
        else:
            v_['KAPPA2'] = args[2]
        if nargin<4:
            v_['KAPPA3'] = 1.0;  #  SIF file default value
        else:
            v_['KAPPA3'] = args[3]
        if nargin<5:
            v_['LB'] = int(5);  #  SIF file default value
        else:
            v_['LB'] = int(args[4])
        if nargin<6:
            v_['UB'] = int(1);  #  SIF file default value
        else:
            v_['UB'] = int(args[5])
        v_['1'] = 1
        v_['MLB'] = -1*v_['LB']
        v_['MUB'] = -1*v_['UB']
        v_['LB+1'] = 1+v_['LB']
        v_['N-UB'] = v_['N']+v_['MUB']
        v_['N-UB-1'] = -1+v_['N-UB']
        v_['-KAPPA3'] = -1.0*v_['KAPPA3']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['LB'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = v_['KAPPA1']+pbm.A[ig,iv]
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
        for I in range(int(v_['LB+1']),int(v_['N-UB-1'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = v_['KAPPA1']+pbm.A[ig,iv]
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
        for I in range(int(v_['N-UB']),int(v_['N'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = v_['KAPPA1']+pbm.A[ig,iv]
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = v_['-KAPPA3']+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'SQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2x_ii( 'CB', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'SQ')
            ielftype = arrset(ielftype, ie, iet_["SQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Q'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'CB')
            ielftype = arrset(ielftype, ie, iet_["CB"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('L2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'L2')
        for I in range(int(v_['1']),int(v_['LB'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['KAPPA2'])
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
        for I in range(int(v_['LB+1']),int(v_['N-UB-1'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+UB'] = I+v_['UB']
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['KAPPA2'])
            for J in range(int(v_['I+1']),int(v_['I+UB'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
        for I in range(int(v_['N-UB']),int(v_['N'])+1):
            v_['I-LB'] = I+v_['MLB']
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['I-LB']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['KAPPA2'])
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,v_['-KAPPA3'])
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def SQ(pbm,nargout,*args):

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
    def CB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def L2(pbm,nargout,*args):

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

