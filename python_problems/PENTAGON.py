from s2mpjlib import *
class  PENTAGON(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PENTAGON
#    *********
# 
#    An approximation to the problem of finding 3 points in a 2D
#    pentagon whose minimal distance is maximal.
# 
#    Source:
#    M.J.D. Powell,
#    " TOLMIN: a Fortran package for linearly constrained
#    optimization problems",
#    Report DAMTP 1989/NA2, University of Cambridge, UK, 1989.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "OLR2-AY-6-15"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PENTAGON'

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
        v_['0'] = 0
        v_['1'] = 1
        v_['3'] = 3
        v_['4'] = 4
        v_['2PI/5'] = 1.2566371
        for J in range(int(v_['0']),int(v_['4'])+1):
            v_['RJ'] = float(J)
            v_['TJ'] = v_['2PI/5']*v_['RJ']
            v_['C'+str(J)] = np.cos(v_['TJ'])
            v_['S'+str(J)] = np.sin(v_['TJ'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['3'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['3'])+1):
            for J in range(int(v_['0']),int(v_['4'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(v_['C'+str(J)])+pbm.A[ig,iv]
                iv = ix_['Y'+str(I)]
                pbm.A[ig,iv] = float(v_['S'+str(J)])+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),1.0)
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(-1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(-1.0)))
        if('Y1' in ix_):
            pb.x0[ix_['Y1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y1']),float(0.0)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(0.0)))
        if('Y2' in ix_):
            pb.x0[ix_['Y2']] = float(-1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y2']),float(-1.0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(1.0)))
        if('Y3' in ix_):
            pb.x0[ix_['Y3']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y3']),float(1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eIDIST', iet_)
        elftv = loaset(elftv,it,0,'XA')
        elftv = loaset(elftv,it,1,'YA')
        elftv = loaset(elftv,it,2,'XB')
        elftv = loaset(elftv,it,3,'YB')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'D12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype, ie, iet_["eIDIST"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'D13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype, ie, iet_["eIDIST"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'D32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype, ie, iet_["eIDIST"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='YB')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D12'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D13'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D32'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
#  LO SOLTN              1.36521631D-04
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-AY-6-15"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eIDIST(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        D = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        D9 = D**9
        D10 = D9*D
        f_   = 1.0/D**8
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -16.0*IV_[0]/D9
            g_[1] = -16.0*IV_[1]/D9
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 16.0*(18.0*IV_[0]*IV_[0]-D)/D10
                H_[0,1] = 288.0*IV_[0]*IV_[1]/D10
                H_[1,0] = H_[0,1]
                H_[1,1] = 16.0*(18.0*IV_[1]*IV_[1]-D)/D10
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

