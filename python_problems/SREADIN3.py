from s2mpjlib import *
class  SREADIN3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SREADIN3
#    *********
# 
#    A nonlinear optimal control problem from Nancy Nichols
#    with a periodic boundary condition and special scaling.
#    This problem arises in tide modelling.
# 
#    SIF input: Nick Gould, July 1991.
# 
#    classification = "OOR2-MN-V-V"
# 
#    Number of discretized points in [0,1]
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER n=6, m=3
# IE N                   5              $-PARAMETER n=12, m=6
# IE N                   50             $-PARAMETER n=102, m=51
# IE N                   100            $-PARAMETER n=202, m=101   original value
# IE N                   500            $-PARAMETER n=1002, m=501
# IE N                   1000           $-PARAMETER n=2002, m=1001
# IE N                   2000           $-PARAMETER n=4002, m=2001
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SREADIN3'

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
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   5000           $-PARAMETER n=10002, m=5001
        v_['PI'] = 3.1415926535
        v_['2PI'] = 2.0*v_['PI']
        v_['A'] = 0.07716
        v_['1/A'] = 1.0/v_['A']
        v_['1/2A'] = 0.5*v_['1/A']
        v_['2A'] = 2.0*v_['A']
        v_['-2A'] = -1.0*v_['2A']
        v_['-1/2A'] = 1.0/v_['-2A']
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = 1.0/v_['RN']
        v_['2/H'] = 2.0*v_['RN']
        v_['H/2'] = 0.5*v_['H']
        v_['1/H'] = 1.0*v_['RN']
        v_['-1/H'] = -1.0*v_['RN']
        v_['0'] = 0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
            xscale = arrset(xscale,iv,v_['RN'])
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('I'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            v_['2PITI'] = v_['2PI']*v_['TI']
            v_['CTI'] = np.cos(v_['2PITI'])
            v_['CCTI'] = v_['CTI']*v_['-1/2A']
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['TI-1'] = v_['RI-1']*v_['H']
            v_['2PITI-1'] = v_['2PI']*v_['TI-1']
            v_['CTI-1'] = np.cos(v_['2PITI-1'])
            v_['CCTI-1'] = v_['CTI-1']*v_['-1/2A']
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['CCTI'])+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['CCTI-1'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('PERIOD',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'PERIOD')
        iv = ix_['X'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['H']))
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
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xlower[ix_['X'+str(I)]] = -0.5
            pb.xupper[ix_['X'+str(I)]] = 0.5
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xlower[ix_['U'+str(I)]] = 0.0
            pb.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.25))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'eENERGY', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'X')
        elftp = loaset(elftp,it,0,'T')
        elftp = loaset(elftp,it,1,'HOVER2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['H']
            ename = 'I'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eENERGY')
            ielftype = arrset(ielftype, ie, iet_["eENERGY"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.25)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.25)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TI']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='HOVER2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['H/2']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            ename = 'NC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.25)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.25)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['1/2A']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ig = ig_['I'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['I'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['I'+str(int(v_['I-1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ig = ig_['C'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['NC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['NC'+str(int(v_['I-1']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        #%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        lxs = len(xscale);
        for j in np.arange(0,min(sA2,pb.n,len(xscale))):
            if not xscale[j] is None and xscale[j] != 0.0 and xscale[j] != 1.0:
                for i in find(pbm.A[:,j],lambda x:x!=0):
                      pbm.A[i,j] = pbm.A[i,j]/xscale[j]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EV_[1]
            g_[1] = pbm.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = pbm.elpar[iel_][0]
                H_[0,1] = H_[1,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eENERGY(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C2PIT = np.cos(2.0*3.141592653589*pbm.elpar[iel_][0])
        f_   = pbm.elpar[iel_][1]*EV_[0]*(EV_[1]-C2PIT)**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][1]*(EV_[1]-C2PIT)**2
            g_[1] = pbm.elpar[iel_][1]*2.0*EV_[0]*(EV_[1]-C2PIT)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = pbm.elpar[iel_][1]*2.0*(EV_[1]-C2PIT)
                H_[0,1] = H_[1,0]
                H_[1,1] = pbm.elpar[iel_][1]*2.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

