from s2mpjlib import *
class  MODBEALENE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MODBEALENE
#    *********
#    A variation on Beale's problem in 2 variables
#    This is a nonlinear equation variant of MODBEALE
# 
#    Source: An adaptation by Ph. Toint of Problem 5 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#89.
#    SIF input: Ph. Toint, Mar 2003.
#               Nick Gould (nonlinear equation version), Jan 2019
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#    The number of variables is  2 * N/2
# 
#           Alternative values for the SIF file parameters:
# IE N/2                 1              $-PARAMETER     original value
# IE N/2                 2              $-PARAMETER
# IE N/2                 5              $-PARAMETER
# IE N/2                 100            $-PARAMETER
# IE N/2                 1000           $-PARAMETER
# IE N/2                 10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MODBEALENE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N/2'] = int(5);  #  SIF file default value
        else:
            v_['N/2'] = int(args[0])
        if nargin<2:
            v_['ALPHA'] = float(50.0);  #  SIF file default value
        else:
            v_['ALPHA'] = float(args[1])
        v_['1'] = 1
        v_['N'] = v_['N/2']+v_['N/2']
        v_['N/2-1'] = -1+v_['N/2']
        v_['ALPHINV'] = 1.0/v_['ALPHA']
        v_['RALPHINV'] = np.sqrt(v_['ALPHINV'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(J),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N/2-1'])+1):
            v_['I-1'] = -1+I
            v_['2I-1'] = v_['I-1']+v_['I-1']
            v_['J'] = 1+v_['2I-1']
            v_['J+1'] = 1+v_['J']
            v_['J+2'] = 2+v_['J']
            [ig,ig_,_] = s2mpj_ii('BA'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'BA'+str(I))
            [ig,ig_,_] = s2mpj_ii('BB'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'BB'+str(I))
            [ig,ig_,_] = s2mpj_ii('BC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'BC'+str(I))
            [ig,ig_,_] = s2mpj_ii('L'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'L'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['J+1']))]])
            valA = np.append(valA,float(6.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['J+2']))]])
            valA = np.append(valA,float(-1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['RALPHINV']))
        [ig,ig_,_] = s2mpj_ii('BA'+str(int(v_['N/2'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'BA'+str(int(v_['N/2'])))
        [ig,ig_,_] = s2mpj_ii('BB'+str(int(v_['N/2'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'BB'+str(int(v_['N/2'])))
        [ig,ig_,_] = s2mpj_ii('BC'+str(int(v_['N/2'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'BC'+str(int(v_['N/2'])))
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
        for I in range(int(v_['1']),int(v_['N/2'])+1):
            self.gconst = arrset(self.gconst,ig_['BA'+str(I)],float(1.5))
            self.gconst = arrset(self.gconst,ig_['BB'+str(I)],float(2.25))
            self.gconst = arrset(self.gconst,ig_['BC'+str(I)],float(2.625))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePRODB', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'POW')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N/2'])+1):
            v_['I-1'] = -1+I
            v_['2I-1'] = v_['I-1']+v_['I-1']
            v_['J'] = 1+v_['2I-1']
            v_['J+1'] = 1+v_['J']
            ename = 'AE'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'ePRODB')
                ielftype = arrset(ielftype,ie,iet_['ePRODB'])
            vname = 'X'+str(int(v_['J']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['J+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='POW')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            ename = 'BE'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'ePRODB')
                ielftype = arrset(ielftype,ie,iet_['ePRODB'])
            vname = 'X'+str(int(v_['J']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['J+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='POW')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(2.0))
            ename = 'CE'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'ePRODB')
                ielftype = arrset(ielftype,ie,iet_['ePRODB'])
            vname = 'X'+str(int(v_['J']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['J+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='POW')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(3.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N/2'])+1):
            ig = ig_['BA'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['AE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['BB'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['BE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['BC'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN                0.0
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
    def ePRODB(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        T = 1.0-EV_[1]**self.elpar[iel_][0]
        POWM1 = self.elpar[iel_][0]-1.0
        W = -self.elpar[iel_][0]*EV_[1]**POWM1
        f_   = EV_[0]*T
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = T
            g_[1] = EV_[0]*W
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 0.0
                H_[0,1] = W
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      -EV_[0]*self.elpar[iel_][0]*POWM1*EV_[1]**(self.elpar[iel_][0]-2.0))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

