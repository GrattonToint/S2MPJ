from s2mpjlib import *
class  MANNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MANNE
#    *********
# 
#    A variable dimension econometric equilibrium problem
#    suggested by A. Manne
# 
#    Source:
#    B. Murtagh and M. Saunders,
#    Mathematical Programming Studies 16, pp. 84-117,
#    (example 5.12).
# 
#    SIF input: N. Gould and Ph. Toint, March 1990.
# 
#    classification = "OOR2-MN-V-V"
# 
#    Number of periods
#    The number of variables in the problem N = 3*T
# 
#           Alternative values for the SIF file parameters:
# IE T                   100            $-PARAMETER n = 300    original value
# IE T                   365            $-PARAMETER n = 995
# IE T                   1000           $-PARAMETER n = 3000
# IE T                   2000           $-PARAMETER n = 6000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MANNE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MANNE'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['T'] = int(4);  #  SIF file default value
        else:
            v_['T'] = int(args[0])
        v_['GROW'] = 0.03
        v_['BETA'] = 0.95
        v_['XK0'] = 3.0
        v_['XC0'] = 0.95
        v_['XI0'] = 0.05
        v_['B'] = 0.25
        v_['BPROB'] = 1.0
        v_['1'] = 1
        v_['2'] = 2
        v_['T-1'] = -1+v_['T']
        v_['T-2'] = -2+v_['T']
        v_['LOGXK'] = np.log(v_['XK0'])
        v_['BLOGX'] = v_['LOGXK']*v_['B']
        v_['XK0**B'] = np.exp(v_['BLOGX'])
        v_['NUM'] = v_['XC0']+v_['XI0']
        v_['A'] = v_['NUM']/v_['XK0**B']
        v_['1-B'] = 1.0-v_['B']
        v_['1+G'] = 1.0+v_['GROW']
        v_['LOG1+G'] = np.log(v_['1+G'])
        v_['SOME'] = v_['LOG1+G']*v_['1-B']
        v_['GFAC'] = np.exp(v_['SOME'])
        v_['AT1'] = v_['A']*v_['GFAC']
        v_['BT1'] = 0.0+v_['BETA']
        for J in range(int(v_['2']),int(v_['T'])+1):
            v_['J-1'] = -1+J
            v_['AT'+str(J)] = v_['AT'+str(int(v_['J-1']))]*v_['GFAC']
            v_['BT'+str(J)] = v_['BT'+str(int(v_['J-1']))]*v_['BETA']
        v_['1-BETA'] = 1.0-v_['BETA']
        v_['1/1-BETA'] = 1.0/v_['1-BETA']
        v_['BT'+str(int(v_['T']))] = v_['BT'+str(int(v_['T']))]*v_['1/1-BETA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'C'+str(I))
            [iv,ix_,_] = s2mpj_ii('I'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'I'+str(I))
            [iv,ix_,_] = s2mpj_ii('K'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'K'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['T'])+1):
            [ig,ig_,_] = s2mpj_ii('NL'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'NL'+str(I))
            iv = ix_['C'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['I'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['T-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('L'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'L'+str(I))
            iv = ix_['K'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['K'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['I'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L'+str(int(v_['T'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L'+str(int(v_['T'])))
        iv = ix_['K'+str(int(v_['T']))]
        pbm.A[ig,iv] = float(v_['GROW'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('L'+str(int(v_['T'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L'+str(int(v_['T'])))
        iv = ix_['I'+str(int(v_['T']))]
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
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['K1']] = 3.05
        pb.xupper[ix_['K1']] = 3.05
        for I in range(int(v_['2']),int(v_['T'])+1):
            pb.xlower[ix_['K'+str(I)]] = 3.05
        v_['1.04**T'] = 0.05
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['1.04**T'] = 1.04*v_['1.04**T']
            pb.xlower[ix_['C'+str(I)]] = 0.95
            pb.xlower[ix_['I'+str(I)]] = 0.05
            pb.xupper[ix_['I'+str(I)]] = v_['1.04**T']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('K1' in ix_):
            pb.x0[ix_['K1']] = float(3.05)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['K1']),float(3.05)))
        for I in range(int(v_['2']),int(v_['T'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['I-1/10'] = 0.1*v_['RI-1']
            v_['VAL'] = 3.0+v_['I-1/10']
            if('K'+str(I) in ix_):
                pb.x0[ix_['K'+str(I)]] = float(v_['VAL'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['K'+str(I)]),float(v_['VAL'])))
        for I in range(int(v_['1']),int(v_['T'])+1):
            if('C'+str(I) in ix_):
                pb.x0[ix_['C'+str(I)]] = float(0.95)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['C'+str(I)]),float(0.95)))
            if('I'+str(I) in ix_):
                pb.x0[ix_['I'+str(I)]] = float(0.05)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['I'+str(I)]),float(0.05)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eLOGS', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'ePOWER', iet_)
        elftv = loaset(elftv,it,0,'K')
        elftp = []
        elftp = loaset(elftp,it,0,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['T'])+1):
            ename = 'LOGC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eLOGS')
            ielftype = arrset(ielftype, ie, iet_["eLOGS"])
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='C')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'KS'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
            ielftype = arrset(ielftype, ie, iet_["ePOWER"])
            vname = 'K'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='K')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='B')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['LOGC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['BT'+str(I)]))
            ig = ig_['NL'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['KS'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['AT'+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -9.7457259D-01
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eLOGS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.log(EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/EV_[0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOWER(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**pbm.elpar[iel_][0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EV_[0]**(pbm.elpar[iel_][0]-1.0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0]  = (
                      pbm.elpar[iel_][0]*(pbm.elpar[iel_][0]-1.0)*EV_[0]**(pbm.elpar[iel_][0]-2.0))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

