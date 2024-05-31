from s2xlib import *
class  ODFITS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple Origin/Destination matrix fit using a minimum entropy
#    approach.  The objective is a combination of different aims, namely
#    to be close to an a priori matrix for some entries, to be consistent
#    with some traffic counts (for some entries) and to be small (for entries
#    where nothing else is known).
# 
#    The objective function is of the form
#         SUM   m T [ ln( T / a ) - 1 ] + E   SUM  T [ ln ( T  ) - 1 ]
#        i in I  i i       i   i            i in J  i        i
#                +  g   SUM   q  F [ ln( F / c ) - 1 ]
#                     i in K   i  i       i   i
#    with the constraints that all Ti and Fi be positive and that
#                         F  =  SUM p   T
#                          i     j   ij  j
#    where the pij represent path weights from an a priori assignment.
# 
#    Source: a modification of an example in
#    L.G. Willumsen,
#    "Origin-Destination Matrix: static estimation"
#    in "Concise Encyclopedia of Traffic and Transportation Systems"
#    (M. Papageorgiou, ed.), Pergamon Press, 1991.
# 
#    M. Bierlaire, private communication, 1991.
# 
#    SIF input: Ph Toint, Dec 1991.
# 
#    classification = "OLR2-MN-10-6"
# 
#    Number of available traffic counts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ODFITS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'ODFITS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['ARCS'] = 6
        v_['TC1'] = 100.0
        v_['TC2'] = 500.0
        v_['TC3'] = 400.0
        v_['TC4'] = 1100.0
        v_['TC5'] = 600.0
        v_['TC6'] = 700.0
        v_['QLT1'] = 1.0
        v_['QLT2'] = 1.0
        v_['QLT3'] = 1.0
        v_['QLT4'] = 1.0
        v_['QLT5'] = 1.0
        v_['QLT6'] = 1.0
        v_['P131'] = 1.0
        v_['P132'] = 0.0
        v_['P133'] = 0.0
        v_['P134'] = 0.0
        v_['P135'] = 0.0
        v_['P136'] = 0.0
        v_['P141'] = 0.0
        v_['P142'] = 1.0
        v_['P143'] = 0.0
        v_['P144'] = 1.0
        v_['P145'] = 0.0
        v_['P146'] = 0.0
        v_['P231'] = 0.0
        v_['P232'] = 0.0
        v_['P233'] = 1.0
        v_['P234'] = 1.0
        v_['P235'] = 1.0
        v_['P236'] = 0.0
        v_['P241'] = 0.0
        v_['P242'] = 0.0
        v_['P243'] = 0.0
        v_['P244'] = 1.0
        v_['P245'] = 1.0
        v_['P246'] = 1.0
        v_['APV13'] = 90.0
        v_['APV14'] = 450.0
        v_['APV23'] = 360.0
        v_['MU13'] = 0.5
        v_['MU14'] = 0.5
        v_['MU23'] = 0.5
        v_['1/MU13'] = 1.0/v_['MU13']
        v_['1/MU14'] = 1.0/v_['MU14']
        v_['1/MU23'] = 1.0/v_['MU23']
        v_['GAMMA'] = 1.5
        v_['ENTROP'] = 0.2
        v_['1/ENTR'] = 1.0/v_['ENTROP']
        v_['1'] = 1
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            v_['1/QLT'+str(I)] = 1.0/v_['QLT'+str(I)]
            v_['G/QLT'+str(I)] = v_['1/QLT'+str(I)]*v_['GAMMA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('T13',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T13')
        [iv,ix_,_] = s2x_ii('T14',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T14')
        [iv,ix_,_] = s2x_ii('T23',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T23')
        [iv,ix_,_] = s2x_ii('T24',ix_)
        pb.xnames=arrset(pb.xnames,iv,'T24')
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [iv,ix_,_] = s2x_ii('F'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'F'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('AP13',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['T13']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/MU13']))
        [ig,ig_,_] = s2x_ii('AP14',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['T14']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/MU14']))
        [ig,ig_,_] = s2x_ii('AP23',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['T23']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/MU23']))
        [ig,ig_,_] = s2x_ii('AP24',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['T24']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/ENTR']))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [ig,ig_,_] = s2x_ii('CP'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['F'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['G/QLT'+str(I)]))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [ig,ig_,_] = s2x_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['F'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['T13']
            pbm.A[ig,iv] = float(v_['P13'+str(I)])+pbm.A[ig,iv]
            iv = ix_['T14']
            pbm.A[ig,iv] = float(v_['P14'+str(I)])+pbm.A[ig,iv]
            iv = ix_['T23']
            pbm.A[ig,iv] = float(v_['P23'+str(I)])+pbm.A[ig,iv]
            iv = ix_['T24']
            pbm.A[ig,iv] = float(v_['P24'+str(I)])+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.1)
        pb.xupper = np.full((pb.n,1),+float('inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['T13']] = float(v_['APV13'])
        pb.x0[ix_['T14']] = float(v_['APV14'])
        pb.x0[ix_['T23']] = float(v_['APV23'])
        pb.x0[ix_['T24']] = float(1.0)
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            pb.x0[ix_['F'+str(I)]] = float(v_['TC'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eXLOGX', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'DEN')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'TFIT13'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
        vname = 'T13'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.1,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='DEN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['APV13']))
        ename = 'TFIT23'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
        vname = 'T23'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.1,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='DEN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['APV23']))
        ename = 'TFIT14'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
        vname = 'T14'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.1,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='DEN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['APV14']))
        ename = 'TFIT24'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
        vname = 'T24'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.1,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='DEN')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            ename = 'CFIT'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
            ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
            vname = 'F'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,0.1,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='DEN')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TC'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['AP13']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TFIT13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AP14']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TFIT14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AP23']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TFIT23'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['AP24']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TFIT24'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            ig = ig_['CP'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CFIT'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
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
        pb.pbclass = "OLR2-MN-10-6"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXLOGX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGX = np.log(EV_[0]/pbm.elpar[iel_][0])
        f_   = EV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0+LOGX
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0/EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

