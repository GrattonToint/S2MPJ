from s2xlib import *
class  HALDMADS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HALDMADS
#    *********
# 
#    A nonlinear minmax problem in five variables.
# 
#    Source: 
#    J. Hald and K. Madsen,
#    "Combined LP and quasi-Newton methods for minimax optimization",
#    Mathematical Programming 20, pp. 49-62, 1981.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "LOR2-AN-6-42"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HALDMADS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HALDMADS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['5'] = 5
        v_['21'] = 21
        v_['T'] = -1.0
        for I in range(int(v_['1']),int(v_['21'])+1):
            v_['Y'+str(I)] = v_['T']
            v_['EY'+str(I)] = np.exp(v_['T'])
            v_['T'] = 0.1+v_['T']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['5'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2x_ii('U',ix_)
        pb.xnames=arrset(pb.xnames,iv,'U')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['21'])+1):
            [ig,ig_,_] = s2x_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(I))
            iv = ix_['U']
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('MF'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'MF'+str(I))
            iv = ix_['U']
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['21'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],float(v_['EY'+str(I)]))
            v_['-EY'+str(I)] = -1.0*v_['EY'+str(I)]
            pbm.gconst = arrset(pbm.gconst,ig_['MF'+str(I)],float(v_['-EY'+str(I)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(0.5)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eHM', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftp = []
        elftp = loaset(elftp,it,0,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['21'])+1):
            ename = 'EL'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eHM')
            ielftype = arrset(ielftype, ie, iet_["eHM"])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='Y')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['21'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            ig = ig_['MF'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "LOR2-AN-6-42"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHM(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        YY = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]
        YYY = YY*pbm.elpar[iel_][0]
        N = EV_[0]+pbm.elpar[iel_][0]*EV_[1]
        D = 1.0+EV_[2]*pbm.elpar[iel_][0]+EV_[3]*YY+EV_[4]*YYY
        DD = D*D
        DDD = DD*D
        f_   = N/D
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/D
            g_[1] = pbm.elpar[iel_][0]/D
            g_[2] = -N*pbm.elpar[iel_][0]/DD
            g_[3] = -N*YY/DD
            g_[4] = -N*YYY/DD
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,2] = -pbm.elpar[iel_][0]/DD
                H_[2,0] = H_[0,2]
                H_[0,3] = -YY/DD
                H_[3,0] = H_[0,3]
                H_[0,4] = -YYY/DD
                H_[4,0] = H_[0,4]
                H_[1,2] = -YY/DD
                H_[2,1] = H_[1,2]
                H_[1,3] = -YYY/DD
                H_[3,1] = H_[1,3]
                H_[1,4] = -YY*YY/DD
                H_[4,1] = H_[1,4]
                H_[2,2] = 2.0*N*YY/DDD
                H_[2,3] = 2.0*N*YYY/DDD
                H_[3,2] = H_[2,3]
                H_[2,4] = 2.0*N*YY*YY/DDD
                H_[4,2] = H_[2,4]
                H_[3,3] = 2.0*N*YY*YY/DDD
                H_[3,4] = 2.0*N*YY*YYY/DDD
                H_[4,3] = H_[3,4]
                H_[4,4] = 2.0*N*YYY*YYY/DDD
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

