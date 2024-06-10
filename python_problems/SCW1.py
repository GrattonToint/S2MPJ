from s2mpjlib import *
class  SCW1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SCW1
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
#    posits a piecewise constant solution that oscilates between
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

    name = 'SCW1'

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
            v_['K'] = int(7);  #  SIF file default value
        else:
            v_['K'] = int(args[0])
# IE K                   3              $-PARAMETER
# IE K                   7              $-PARAMETER     original value
# IE K                   10             $-PARAMETER
# IE K                   100            $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['K+1'] = 1+v_['K']
        v_['RK+1'] = float(v_['K+1'])
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['2PI'] = 2.0*v_['PI']
        v_['2PI/K+1'] = v_['2PI']/v_['RK+1']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['K+1'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('S',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['K'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('CON'+str(I),ig_)
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['T'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['T'+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['K'])+1):
            pb.xlower[ix_['T'+str(I)]] = 0.0
            pb.xupper[ix_['T'+str(I)]] = v_['2PI']
        pb.xlower[ix_['T'+str(int(v_['K+1']))]] = v_['2PI']
        pb.xupper[ix_['T'+str(int(v_['K+1']))]] = v_['2PI']
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
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSINT', iet_)
        elftv = loaset(elftv,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'eCOST', iet_)
        elftv = loaset(elftv,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['K+1'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSINT')
            ielftype = arrset(ielftype, ie, iet_["eSINT"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOST')
            ielftype = arrset(ielftype, ie, iet_["eCOST"])
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gMAXSQ',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['C']
        pbm.grftype = arrset(pbm.grftype,ig,'gMAXSQ')
        for I in range(int(v_['0']),int(v_['K'])+1,int(v_['2'])):
            v_['I+1'] = 1+I
            ig = ig_['C']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        for I in range(int(v_['1']),int(v_['K'])+1,int(v_['2'])):
            v_['I+1'] = 1+I
            ig = ig_['C']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        ig = ig_['S']
        pbm.grftype = arrset(pbm.grftype,ig,'gMAXSQ')
        for I in range(int(v_['0']),int(v_['K'])+1,int(v_['2'])):
            v_['I+1'] = 1+I
            ig = ig_['S']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        for I in range(int(v_['1']),int(v_['K'])+1,int(v_['2'])):
            v_['I+1'] = 1+I
            ig = ig_['S']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SCW                 0.0
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSINT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S = np.sin(EV_[0])
        f_   = S
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
                H_[0,0] = -S
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOST(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C = np.cos(EV_[0])
        f_   = C
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -np.sin(EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -C
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

