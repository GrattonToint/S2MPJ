from s2xlib import *
class  SCW2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SCW2
#    *********
# 
#    Source: a discretization of an infinite-demsional problem proposed 
#    by Simon Chandler-Wilde (U. Reading):
# 
#    Given a function u in C[0,2 pi] with ||u||_infty <= 1, find the 
#    supremum of c^2(u) + s^2(u), where
#      c(u) = int_0^2 pi cos(t)u(t) dt and
#      s(u) = int_0^2 pi sin(t)u(t) dt      
# 
#    The discretized version ignores the required continuity, and 
#    posits a piecewise constant solution that varies anywhere between
#    plus and minus one. The anticipated solution is -16.
# 
#    SIF input: Nick Gould, July 2020
# 
#    classification = "SLR2-MN-V-V"
# 
#    Number of internal knots
# 
#           Alternative values for the SIF file parameters:
# IE K                   1              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SCW2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'SCW2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['K'] = int(10);  #  SIF file default value
        else:
            v_['K'] = int(args[0])
# IE K                   100            $-PARAMETER
# IE K                   1000           $-PARAMETER     original value
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['ONE'] = 1.0
        v_['K+1'] = 1+v_['K']
        v_['RK'] = float(v_['K'])
        v_['RK+1'] = float(v_['K+1'])
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['2PI'] = 2.0*v_['PI']
        v_['2PI/K+1'] = v_['2PI']/v_['RK+1']
        v_['1/K'] = v_['ONE']/v_['RK']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['K+1'])+1):
            [iv,ix_,_] = s2x_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
        for I in range(int(v_['0']),int(v_['K'])+1):
            [iv,ix_,_] = s2x_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('S',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('C',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['K'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('CON'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CON'+str(I))
            iv = ix_['T'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['T'+str(I)]
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['T'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['T'+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['K'])+1):
            pb.xlower[ix_['T'+str(I)]] = 0.0
            pb.xupper[ix_['T'+str(I)]] = v_['2PI']
        pb.xlower[ix_['T'+str(int(v_['K+1']))]] = v_['2PI']
        pb.xupper[ix_['T'+str(int(v_['K+1']))]] = v_['2PI']
        for I in range(int(v_['0']),int(v_['K'])+1):
            pb.xlower[ix_['U'+str(I)]] = -1.0
            pb.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('T'+str(int(v_['1'])) in ix_):
            pb.x0[ix_['T'+str(int(v_['1']))]] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['T'+str(int(v_['1']))]),float(1.0)))
        for I in range(int(v_['2']),int(v_['K'])+1):
            v_['RI'] = float(I)
            v_['START'] = v_['RI']*v_['2PI/K+1']
            if('T'+str(I) in ix_):
                pb.x0[ix_['T'+str(I)]] = float(0.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['T'+str(I)]),float(0.0)))
        for I in range(int(v_['0']),int(v_['K'])+1):
            v_['RI'] = float(I)
            v_['START'] = v_['RI']*v_['1/K']
            if('U'+str(I) in ix_):
                pb.x0[ix_['U'+str(I)]] = float(v_['START'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['U'+str(I)]),float(v_['START'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eUSINT', iet_)
        elftv = loaset(elftv,it,0,'T')
        elftv = loaset(elftv,it,1,'U')
        [it,iet_,_] = s2x_ii( 'eUCOST', iet_)
        elftv = loaset(elftv,it,0,'T')
        elftv = loaset(elftv,it,1,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['K'])+1):
            v_['I+1'] = 1+I
            ename = 'US'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eUSINT')
            ielftype = arrset(ielftype, ie, iet_["eUSINT"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'USP'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eUSINT')
            ielftype = arrset(ielftype, ie, iet_["eUSINT"])
            vname = 'T'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'UC'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eUCOST')
            ielftype = arrset(ielftype, ie, iet_["eUCOST"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'UCP'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eUCOST')
            ielftype = arrset(ielftype, ie, iet_["eUCOST"])
            vname = 'T'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gMAXSQ',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['C']
        pbm.grftype = arrset(pbm.grftype,ig,'gMAXSQ')
        for I in range(int(v_['0']),int(v_['K'])+1):
            ig = ig_['C']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['USP'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['US'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['S']
        pbm.grftype = arrset(pbm.grftype,ig,'gMAXSQ')
        for I in range(int(v_['0']),int(v_['K'])+1):
            ig = ig_['S']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['UCP'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['UC'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "SLR2-MN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eUSINT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S = np.sin(EV_[0])
        C = np.cos(EV_[0])
        f_   = EV_[1]*S
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*C
            g_[1] = S
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = -EV_[1]*S
                H_[0,1] = C
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eUCOST(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S = np.sin(EV_[0])
        C = np.cos(EV_[0])
        f_   = EV_[1]*C
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EV_[1]*S
            g_[1] = C
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = -EV_[1]*C
                H_[0,1] = -S
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gMAXSQ(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= -GVAR_*GVAR_
        if nargout>1:
            g_ = -GVAR_-GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

