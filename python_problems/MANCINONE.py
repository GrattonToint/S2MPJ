from s2mpjlib import *
class  MANCINONE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MANCINONE
#    *********
# 
#    Mancino's function with variable dimension.
#    This is a nonlinear equation variant of MANCINO
# 
#    Source:
#    E. Spedicato,
#    "Computational experience with quasi-Newton algorithms for
#    minimization problems of moderate size",
#    Report N-175, CISE, Milano, 1975.
# 
#    See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995.
#               Nick Gould (nonlinear equation version), Jan 2019
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#    The definitions
#      s_{i,j} = \sin \log v_{i,j}   and s_{i,j} = \cos \log v_{i,j}
#    have been used.  It seems that the additional exponent ALPHA
#    in Buckley is a typo.
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   20             $-PARAMETER
# IE N                   30             $-PARAMETER Schittkowski #391
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MANCINONE'

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
        if nargin<2:
            v_['ALPHA'] = int(5);  #  SIF file default value
        else:
            v_['ALPHA'] = int(args[1])
        if nargin<3:
            v_['BETA'] = float(14.0);  #  SIF file default value
        else:
            v_['BETA'] = float(args[2])
        if nargin<4:
            v_['GAMMA'] = int(3);  #  SIF file default value
        else:
            v_['GAMMA'] = int(args[3])
        v_['RALPHA'] = float(v_['ALPHA'])
        v_['RN'] = float(v_['N'])
        v_['N-1'] = -1+v_['N']
        v_['RN-1'] = float(v_['N-1'])
        v_['N-1SQ'] = v_['RN-1']*v_['RN-1']
        v_['BETAN'] = v_['BETA']*v_['RN']
        v_['BETAN2'] = v_['BETAN']*v_['BETAN']
        v_['AL+1'] = 1.0+v_['RALPHA']
        v_['A1SQ'] = v_['AL+1']*v_['AL+1']
        v_['F0'] = v_['A1SQ']*v_['N-1SQ']
        v_['F1'] = -1.0*v_['F0']
        v_['F2'] = v_['BETAN2']+v_['F1']
        v_['F3'] = 1.0/v_['F2']
        v_['F4'] = v_['BETAN']*v_['F3']
        v_['A'] = -1.0*v_['F4']
        v_['-N/2'] = -0.5*v_['RN']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(v_['BETAN']))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['I-N/2'] = v_['RI']+v_['-N/2']
            v_['CI'] = 1.0
            for J in range(int(v_['1']),int(v_['GAMMA'])+1):
                v_['CI'] = v_['CI']*v_['I-N/2']
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['CI']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['H'] = 0.0
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RJ'] = float(J)
                v_['1/J'] = 1.0/v_['RJ']
                v_['I/J'] = v_['RI']*v_['1/J']
                v_['SQI/J'] = np.sqrt(v_['I/J'])
                v_['LIJ'] = np.log(v_['SQI/J'])
                v_['SIJ'] = np.sin(v_['LIJ'])
                v_['CIJ'] = np.cos(v_['LIJ'])
                v_['SA'] = 1.0
                v_['CA'] = 1.0
                for K in range(int(v_['1']),int(v_['ALPHA'])+1):
                    v_['SA'] = v_['SA']*v_['SIJ']
                    v_['CA'] = v_['CA']*v_['CIJ']
                v_['SCA'] = v_['SA']+v_['CA']
                v_['HIJ'] = v_['SQI/J']*v_['SCA']
                v_['H'] = v_['H']+v_['HIJ']
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['1/J'] = 1.0/v_['RJ']
                v_['I/J'] = v_['RI']*v_['1/J']
                v_['SQI/J'] = np.sqrt(v_['I/J'])
                v_['LIJ'] = np.log(v_['SQI/J'])
                v_['SIJ'] = np.sin(v_['LIJ'])
                v_['CIJ'] = np.cos(v_['LIJ'])
                v_['SA'] = 1.0
                v_['CA'] = 1.0
                for K in range(int(v_['1']),int(v_['ALPHA'])+1):
                    v_['SA'] = v_['SA']*v_['SIJ']
                    v_['CA'] = v_['CA']*v_['CIJ']
                v_['SCA'] = v_['SA']+v_['CA']
                v_['HIJ'] = v_['SQI/J']*v_['SCA']
                v_['H'] = v_['H']+v_['HIJ']
            v_['I-N/2'] = v_['RI']+v_['-N/2']
            v_['CI'] = 1.0
            for J in range(int(v_['1']),int(v_['GAMMA'])+1):
                v_['CI'] = v_['CI']*v_['I-N/2']
            v_['TMP'] = v_['H']+v_['CI']
            v_['XI0'] = v_['TMP']*v_['A']
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(v_['XI0'])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['XI0'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eMANC', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'II')
        elftp = loaset(elftp,it,1,'JJ')
        elftp = loaset(elftp,it,2,'AL')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RJ'] = float(J)
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eMANC')
                ielftype = arrset(ielftype,ie,iet_["eMANC"])
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='II')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RI']))
                posep = np.where(elftp[ielftype[ie]]=='JJ')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RJ']))
                posep = np.where(elftp[ielftype[ie]]=='AL')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RALPHA']))
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eMANC')
                ielftype = arrset(ielftype,ie,iet_["eMANC"])
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='II')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RI']))
                posep = np.where(elftp[ielftype[ie]]=='JJ')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RJ']))
                posep = np.where(elftp[ielftype[ie]]=='AL')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['RALPHA']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['G'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eMANC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        IAL = self.elpar[iel_][2]
        IA1 = IAL-1
        A2 = self.elpar[iel_][2]-2.0
        IA2 = IAL-2
        IA3 = IAL-3
        INVIJ = EV_[0,0]*EV_[0,0]+self.elpar[iel_][0]/self.elpar[iel_][1]
        VIJ = np.sqrt(INVIJ)
        V2 = VIJ*VIJ
        DVIJ = EV_[0,0]/VIJ
        LIJ = np.log(VIJ)
        SIJ = np.sin(LIJ)
        CIJ = np.cos(LIJ)
        DSDX = CIJ*DVIJ/VIJ
        DCDX = -SIJ*DVIJ/VIJ
        SUMAL = SIJ**IAL+CIJ**IAL
        DSUMAL = self.elpar[iel_][2]*(DSDX*SIJ**IA1+DCDX*CIJ**IA1)
        SCIJ = SIJ*CIJ
        DSCIJ = SIJ*DCDX+DSDX*CIJ
        SAL = SIJ**IA2-CIJ**IA2
        DSAL = A2*(DSDX*SIJ**IA3-DCDX*CIJ**IA3)
        B = SUMAL+self.elpar[iel_][2]*SCIJ*SAL
        DBDX = DSUMAL+self.elpar[iel_][2]*(DSCIJ*SAL+SCIJ*DSAL)
        f_   = VIJ*SUMAL
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]*B/VIJ
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = (B+EV_[0,0]*DBDX)/VIJ-B*EV_[0,0]*DVIJ/V2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

