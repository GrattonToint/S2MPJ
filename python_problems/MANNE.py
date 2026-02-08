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
#    classification = "C-COOR2-MN-V-V"
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
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MANNE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['T'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'C'+str(I))
            [iv,ix_,_] = s2mpj_ii('I'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'I'+str(I))
            [iv,ix_,_] = s2mpj_ii('K'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'K'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['T'])+1):
            [ig,ig_,_] = s2mpj_ii('NL'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'NL'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['C'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['I'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['T-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('L'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'L'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['K'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['K'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['I'+str(I)]])
            valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('L'+str(int(v_['T'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L'+str(int(v_['T'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['K'+str(int(v_['T']))]])
        valA = np.append(valA,float(v_['GROW']))
        [ig,ig_,_] = s2mpj_ii('L'+str(int(v_['T'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'L'+str(int(v_['T'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['I'+str(int(v_['T']))]])
        valA = np.append(valA,float(-1.0))
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
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = -np.full((self.nle,1),float('inf'))
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['K1']] = 3.05
        self.xupper[ix_['K1']] = 3.05
        for I in range(int(v_['2']),int(v_['T'])+1):
            self.xlower[ix_['K'+str(I)]] = 3.05
        v_['1.04**T'] = 0.05
        for I in range(int(v_['1']),int(v_['T'])+1):
            v_['1.04**T'] = 1.04*v_['1.04**T']
            self.xlower[ix_['C'+str(I)]] = 0.95
            self.xlower[ix_['I'+str(I)]] = 0.05
            self.xupper[ix_['I'+str(I)]] = v_['1.04**T']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('K1' in ix_):
            self.x0[ix_['K1']] = float(3.05)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['K1']),float(3.05)))
        for I in range(int(v_['2']),int(v_['T'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['I-1/10'] = 0.1*v_['RI-1']
            v_['VAL'] = 3.0+v_['I-1/10']
            if('K'+str(I) in ix_):
                self.x0[ix_['K'+str(I)]] = float(v_['VAL'])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['K'+str(I)]),float(v_['VAL'])))
        for I in range(int(v_['1']),int(v_['T'])+1):
            if('C'+str(I) in ix_):
                self.x0[ix_['C'+str(I)]] = float(0.95)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['C'+str(I)]),float(0.95)))
            if('I'+str(I) in ix_):
                self.x0[ix_['I'+str(I)]] = float(0.05)
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['I'+str(I)]),float(0.05)))
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['T'])+1):
            ename = 'LOGC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eLOGS')
            ielftype = arrset(ielftype,ie,iet_["eLOGS"])
            vname = 'C'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='C')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'KS'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePOWER')
            ielftype = arrset(ielftype,ie,iet_["ePOWER"])
            vname = 'K'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='K')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='B')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['T'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['LOGC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['BT'+str(I)]))
            ig = ig_['NL'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['KS'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['AT'+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -9.7457259D-01
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eLOGS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.log(EV_[0,0])
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/EV_[0,0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOWER(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]**self.elpar[iel_][0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*EV_[0,0]**(self.elpar[iel_][0]-1.0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0]  = (
                      self.elpar[iel_][0]*(self.elpar[iel_][0]-1.0)*EV_[0,0]**(self.elpar[iel_][0]-2.0))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

