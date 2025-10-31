from s2mpjlib import *
class  ROSZMAN1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ROSZMAN1
#    *********
# 
#    NIST Data fitting problem ROSZMAN1 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y =  b1 - b2*x - arctan[b3/(x-b4)]/pi + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#   Reference: Roszman, L., NIST (1979).  
#     Quantum Defects for Sulfur I Atom.
# 
#    classification = "C-CNOR2-MN-4-25"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ROSZMAN1'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 25
        v_['N'] = 4
        v_['1'] = 1
        v_['X1'] = -4868.68
        v_['X2'] = -4868.09
        v_['X3'] = -4867.41
        v_['X4'] = -3375.19
        v_['X5'] = -3373.14
        v_['X6'] = -3372.03
        v_['X7'] = -2473.74
        v_['X8'] = -2472.35
        v_['X9'] = -2469.45
        v_['X10'] = -1894.65
        v_['X11'] = -1893.40
        v_['X12'] = -1497.24
        v_['X13'] = -1495.85
        v_['X14'] = -1493.41
        v_['X15'] = -1208.68
        v_['X16'] = -1206.18
        v_['X17'] = -1206.04
        v_['X18'] = -997.92
        v_['X19'] = -996.61
        v_['X20'] = -996.31
        v_['X21'] = -834.94
        v_['X22'] = -834.66
        v_['X23'] = -710.03
        v_['X24'] = -530.16
        v_['X25'] = -464.17
        v_['Y1'] = 0.252429
        v_['Y2'] = 0.252141
        v_['Y3'] = 0.251809
        v_['Y4'] = 0.297989
        v_['Y5'] = 0.296257
        v_['Y6'] = 0.295319
        v_['Y7'] = 0.339603
        v_['Y8'] = 0.337731
        v_['Y9'] = 0.333820
        v_['Y10'] = 0.389510
        v_['Y11'] = 0.386998
        v_['Y12'] = 0.438864
        v_['Y13'] = 0.434887
        v_['Y14'] = 0.427893
        v_['Y15'] = 0.471568
        v_['Y16'] = 0.461699
        v_['Y17'] = 0.461144
        v_['Y18'] = 0.513532
        v_['Y19'] = 0.506641
        v_['Y20'] = 0.505062
        v_['Y21'] = 0.535648
        v_['Y22'] = 0.533726
        v_['Y23'] = 0.568064
        v_['Y24'] = 0.612886
        v_['Y25'] = 0.624169
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B1']])
            valA = np.append(valA,float(1.0))
            v_['-X'] = -1.0*v_['X'+str(I)]
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['B2']])
            valA = np.append(valA,float(v_['-X']))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('B1' in ix_):
            self.x0[ix_['B1']] = float(0.1)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B1']),float(0.1)))
        if('B2' in ix_):
            self.x0[ix_['B2']] = float(-0.00001)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B2']),float(-0.00001)))
        if('B3' in ix_):
            self.x0[ix_['B3']] = float(1000.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B3']),float(1000.0)))
        if('B4' in ix_):
            self.x0[ix_['B4']] = float(-100.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B4']),float(-100.0)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE7', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE7')
            ielftype = arrset(ielftype,ie,iet_["eE7"])
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               
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
        self.pbclass   = "C-CNOR2-MN-4-25"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,4.0*np.arctan(1.0e0))
        return pbm

    @staticmethod
    def eE7(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V12 = EV_[0]*EV_[0]
        V13 = EV_[0]*V12
        V2MX = EV_[1]-self.elpar[iel_][0]
        V2MX2 = V2MX*V2MX
        V2MX3 = V2MX*V2MX2
        R = V12/V2MX2+1.0
        PIR = self.efpar[0]*R
        PIR2 = PIR*R
        f_   = -np.arctan(EV_[0]/V2MX)/self.efpar[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/(PIR*V2MX)
            g_[1] = EV_[0]/(PIR*V2MX2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*EV_[0]/(PIR2*V2MX3)
                H_[0,1] = 1.0/(PIR*V2MX2)-2.0*V12/(PIR2*V2MX**4)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*V13/(PIR2*V2MX**5)-2.0*EV_[0]/(PIR*V2MX3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

