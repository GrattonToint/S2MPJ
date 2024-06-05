from s2mpjlib import *
class  MINSURF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINSURF
#    *********
#    Variable dimension full rank linear problem
#    A version of the minimum surface problem
#    on the unit square with simple boundary conditions.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "OXR2-MY-64-0"
# 
#    Discretization parameter
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MINSURF'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MINSURF'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['P'] = 7
        v_['1'] = 1
        v_['P+1'] = 1+v_['P']
        v_['RP'] = float(v_['P'])
        v_['RPSQ'] = v_['RP']*v_['RP']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for i in range(int(v_['1']),int(v_['P+1'])+1):
            for j in range(int(v_['1']),int(v_['P+1'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(i)+','+str(j),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(i)+','+str(j))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                [ig,ig_,_] = s2mpj_ii('S'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                pbm.gscale = arrset(pbm.gscale,ig,float(v_['RPSQ']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['S'+str(i)+','+str(j)],float(-1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        v_['2'] = 2
        for i in range(int(v_['2']),int(v_['P'])+1):
            for j in range(int(v_['2']),int(v_['P'])+1):
                pb.xlower[ix_['X'+str(i)+','+str(j)]] = -float('Inf')
                pb.xupper[ix_['X'+str(i)+','+str(j)]] = +float('Inf')
        for i in range(int(v_['1']),int(v_['P+1'])+1):
            pb.xlower[ix_['X'+str(int(v_['1']))+','+str(i)]] = 1.0
            pb.xupper[ix_['X'+str(int(v_['1']))+','+str(i)]] = 1.0
            pb.xlower[ix_['X'+str(int(v_['P+1']))+','+str(i)]] = 1.0
            pb.xupper[ix_['X'+str(int(v_['P+1']))+','+str(i)]] = 1.0
            pb.xlower[ix_['X'+str(i)+','+str(int(v_['1']))]] = 1.0
            pb.xupper[ix_['X'+str(i)+','+str(int(v_['1']))]] = 1.0
            pb.xlower[ix_['X'+str(i)+','+str(int(v_['P+1']))]] = 1.0
            pb.xupper[ix_['X'+str(i)+','+str(int(v_['P+1']))]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for i in range(int(v_['1']),int(v_['P'])+1):
            v_['i+1'] = 1+i
            for j in range(int(v_['1']),int(v_['P'])+1):
                v_['j+1'] = 1+j
                ename = 'A'+str(i)+','+str(j)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'X'+str(i)+','+str(j)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['i+1']))+','+str(int(v_['j+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='W')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'B'+str(i)+','+str(j)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'X'+str(i)+','+str(int(v_['j+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['i+1']))+','+str(j)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='W')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQROOT',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        v_['WEIGHT'] = 0.5*v_['RPSQ']
        for i in range(int(v_['1']),int(v_['P'])+1):
            for j in range(int(v_['1']),int(v_['P'])+1):
                ig = ig_['S'+str(i)+','+str(j)]
                pbm.grftype = arrset(pbm.grftype,ig,'gSQROOT')
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(i)+','+str(j)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WEIGHT']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(i)+','+str(j)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WEIGHT']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OXR2-MY-64-0"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQROOT(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        SQRAL = np.sqrt(GVAR_)
        f_= SQRAL
        if nargout>1:
            g_ = 0.5/SQRAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.25/(GVAR_*SQRAL)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

