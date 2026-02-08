from s2mpjlib import *
class  SEMICON1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SEMICON1
#    *********
# 
#    The semiconductor problem by Rheinboldt, using a finite difference
#    approximation.
# 
#    Source: problem 10 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#    N  = Number of discretized point inside the interval [a, b]
#    LN = Index of the last negative discretization point
#         (the interest is in the negative part)
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SEMICON1'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE LN                  9              $-PARAMETER     original value
        if nargin<2:
            v_['LN'] = int(9);  #  SIF file default value
        else:
            v_['LN'] = int(args[1])
# IE N                   50             $-PARAMETER
# IE LN                  45             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE LN                  90             $-PARAMETER
# IE N                   500            $-PARAMETER
# IE LN                  450            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE LN                  900            $-PARAMETER
# IE N                   5000           $-PARAMETER
# IE LN                  4500           $-PARAMETER
        if nargin<3:
            v_['LAMBDA'] = float(1.0);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[2])
        v_['A'] = -0.00009
        v_['B'] = 0.00001
        v_['UA'] = 0.0
        v_['UB'] = 700.0
        v_['CA'] = 1.0e12
        v_['CB'] = 1.0e13
        v_['BETA'] = 40.0
        v_['LN+1'] = 1+v_['LN']
        v_['N+1'] = 1+v_['N']
        v_['-A'] = -1.0*v_['A']
        v_['B-A'] = v_['B']+v_['-A']
        v_['RN+1'] = float(v_['N+1'])
        v_['TMP'] = 1.0/v_['RN+1']
        v_['H'] = v_['B-A']*v_['TMP']
        v_['H2'] = v_['H']*v_['H']
        v_['LB'] = v_['LAMBDA']*v_['BETA']
        v_['H2CA'] = v_['H2']*v_['CA']
        v_['H2CB'] = v_['H2']*v_['CB']
        v_['LH2CA'] = v_['LAMBDA']*v_['H2CA']
        v_['LH2CB'] = v_['LAMBDA']*v_['H2CB']
        v_['LUA'] = v_['LAMBDA']*v_['UA']
        v_['LUB'] = v_['LAMBDA']*v_['UB']
        v_['ULW'] = -5.0+v_['LUA']
        v_['UUP'] = 5.0+v_['LUB']
        v_['-LB'] = -1.0*v_['LB']
        v_['-LUB'] = -1.0*v_['LUB']
        v_['-LH2CB'] = -1.0*v_['LH2CB']
        v_['0'] = 0
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N+1'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(-2.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['LN'])+1):
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['LH2CA']))
        for I in range(int(v_['LN+1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['-LH2CB']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),v_['UUP'])
        self.xlower = np.full((self.n,1),v_['ULW'])
        self.xlower[ix_['U'+str(int(v_['0']))]] = v_['LUA']
        self.xupper[ix_['U'+str(int(v_['0']))]] = v_['LUA']
        self.xlower[ix_['U'+str(int(v_['N+1']))]] = v_['LUB']
        self.xupper[ix_['U'+str(int(v_['N+1']))]] = v_['LUB']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        self.x0[ix_['U'+str(int(v_['0']))]] = float(v_['LUA'])
        self.x0[ix_['U'+str(int(v_['N+1']))]] = float(v_['LUB'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eWE1', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'LAC')
        elftp = loaset(elftp,it,1,'LAB')
        elftp = loaset(elftp,it,2,'LU')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eWE1')
            ielftype = arrset(ielftype,ie,iet_["eWE1"])
            vname = 'U'+str(I)
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(v_['ULW']),float(v_['UUP']),float(0.0)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='LAC')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LH2CA']))
            posep = np.where(elftp[ielftype[ie]]=='LAB')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-LB']))
            posep = np.where(elftp[ielftype[ie]]=='LU')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LUA']))
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eWE1')
            ielftype = arrset(ielftype,ie,iet_["eWE1"])
            vname = 'U'+str(I)
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(v_['ULW']),float(v_['UUP']),float(0.0)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='LAC')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['-LH2CB']))
            posep = np.where(elftp[ielftype[ie]]=='LAB')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LB']))
            posep = np.where(elftp[ielftype[ie]]=='LU')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LUB']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
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
        self.pbclass   = "C-CNOR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eWE1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVAL  = (
              self.elpar[iel_][0]*np.exp(self.elpar[iel_][1]*(EV_[0,0]-self.elpar[iel_][2])))
        f_   = FVAL
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][1]*FVAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = self.elpar[iel_][1]*self.elpar[iel_][1]*FVAL
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

