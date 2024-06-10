from s2mpjlib import *
class  MSQRTBLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MSQRTBLS
#    *********
# 
#    The dense matrix square root problem by Nocedal and Liu (Case 1)
# 
#    This is a least-squares variant of problem MSQRTB.
# 
#    Source:  problem 204 (p. 93) in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-AN-V-V"
# 
#    Dimension of the matrix ( at least 3)
# 
#           Alternative values for the SIF file parameters:
# IE P                   3              $-PARAMETER n = 9     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   10             $-PARAMETER n = 100
# IE P                   23             $-PARAMETER n = 529
# IE P                   32             $-PARAMETER n = 1024
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MSQRTBLS'

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
            v_['P'] = int(5);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   70             $-PARAMETER n = 4900
        v_['N'] = v_['P']*v_['P']
        v_['1'] = 1
        v_['K'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = 1.0+v_['K']
                v_['K2'] = v_['K']*v_['K']
                v_['B'+str(I)+','+str(J)] = np.sin(v_['K2'])
        v_['B3,1'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
                for T in range(int(v_['1']),int(v_['P'])+1):
                    v_['PROD'] = v_['B'+str(I)+','+str(T)]*v_['B'+str(T)+','+str(J)]
                    v_['A'+str(I)+','+str(J)] = v_['A'+str(I)+','+str(J)]+v_['PROD']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['G'+str(I)+','+str(J)],float(v_['A'+str(I)+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        v_['K'] = 0.0
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['K'] = 1.0+v_['K']
                v_['K2'] = v_['K']*v_['K']
                v_['SK2'] = np.sin(v_['K2'])
                v_['-4SK2/5'] = -0.8*v_['SK2']
                v_['XIJ'] = v_['B'+str(I)+','+str(J)]+v_['-4SK2/5']
                pb.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['XIJ'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'XIT')
        elftv = loaset(elftv,it,1,'XTJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                for T in range(int(v_['1']),int(v_['P'])+1):
                    ename = 'E'+str(I)+','+str(J)+','+str(T)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    vname = 'X'+str(I)+','+str(T)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='XIT')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(T)+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='XTJ')
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
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                for T in range(int(v_['1']),int(v_['P'])+1):
                    ig = ig_['G'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)+','+str(T)]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
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

