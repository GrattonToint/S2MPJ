from s2mpjlib import *
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
#    classification = "C-COLR2-MN-10-6"
# 
#    Number of available traffic counts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ODFITS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('T13',ix_)
        self.xnames=arrset(self.xnames,iv,'T13')
        [iv,ix_,_] = s2mpj_ii('T14',ix_)
        self.xnames=arrset(self.xnames,iv,'T14')
        [iv,ix_,_] = s2mpj_ii('T23',ix_)
        self.xnames=arrset(self.xnames,iv,'T23')
        [iv,ix_,_] = s2mpj_ii('T24',ix_)
        self.xnames=arrset(self.xnames,iv,'T24')
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [iv,ix_,_] = s2mpj_ii('F'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'F'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('AP13',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T13']])
        valA = np.append(valA,float(-1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['1/MU13']))
        [ig,ig_,_] = s2mpj_ii('AP14',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T14']])
        valA = np.append(valA,float(-1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['1/MU14']))
        [ig,ig_,_] = s2mpj_ii('AP23',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T23']])
        valA = np.append(valA,float(-1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['1/MU23']))
        [ig,ig_,_] = s2mpj_ii('AP24',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T24']])
        valA = np.append(valA,float(-1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['1/ENTR']))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [ig,ig_,_] = s2mpj_ii('CP'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['F'+str(I)]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['G/QLT'+str(I)]))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['F'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T13']])
            valA = np.append(valA,float(v_['P13'+str(I)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T14']])
            valA = np.append(valA,float(v_['P14'+str(I)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T23']])
            valA = np.append(valA,float(v_['P23'+str(I)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T24']])
            valA = np.append(valA,float(v_['P24'+str(I)]))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.1)
        self.xupper = np.full((self.n,1),+float('inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['T13']] = float(v_['APV13'])
        self.x0[ix_['T14']] = float(v_['APV14'])
        self.x0[ix_['T23']] = float(v_['APV23'])
        self.x0[ix_['T24']] = float(1.0)
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            self.x0[ix_['F'+str(I)]] = float(v_['TC'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXLOGX', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'DEN')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'TFIT13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
        vname = 'T13'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='DEN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['APV13']))
        ename = 'TFIT23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
        vname = 'T23'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='DEN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['APV23']))
        ename = 'TFIT14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
        vname = 'T14'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='DEN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['APV14']))
        ename = 'TFIT24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXLOGX')
        ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
        vname = 'T24'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='DEN')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            ename = 'CFIT'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eXLOGX')
            ielftype = arrset(ielftype,ie,iet_["eXLOGX"])
            vname = 'F'+str(I)
            [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='DEN')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TC'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['AP13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TFIT13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AP14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TFIT14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AP23']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TFIT23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AP24']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TFIT24'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['ARCS'])+1):
            ig = ig_['CP'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CFIT'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO ODFITS             -2380.026775
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COLR2-MN-10-6"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXLOGX(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGX = np.log(EV_[0]/self.elpar[iel_][0])
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

