from s2mpjlib import *
class  HS57(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS57
#    *********
# 
#    Source: problem 57 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, April 1990
# 
#    classification = "C-CSQR2-AN-2-1"
# 
#    Problem parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS57'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['A1'] = 8.0
        v_['A2'] = 8.0
        v_['A3'] = 10.0
        v_['A4'] = 10.0
        v_['A5'] = 10.0
        v_['A6'] = 10.0
        v_['A7'] = 12.0
        v_['A8'] = 12.0
        v_['A9'] = 12.0
        v_['A10'] = 12.0
        v_['A11'] = 14.0
        v_['A12'] = 14.0
        v_['A13'] = 14.0
        v_['A14'] = 16.0
        v_['A15'] = 16.0
        v_['A16'] = 16.0
        v_['A17'] = 18.0
        v_['A18'] = 18.0
        v_['A19'] = 20.0
        v_['A20'] = 20.0
        v_['A21'] = 20.0
        v_['A22'] = 22.0
        v_['A23'] = 22.0
        v_['A24'] = 22.0
        v_['A25'] = 24.0
        v_['A26'] = 24.0
        v_['A27'] = 24.0
        v_['A28'] = 26.0
        v_['A29'] = 26.0
        v_['A30'] = 26.0
        v_['A31'] = 28.0
        v_['A32'] = 28.0
        v_['A33'] = 30.0
        v_['A34'] = 30.0
        v_['A35'] = 30.0
        v_['A36'] = 32.0
        v_['A37'] = 32.0
        v_['A38'] = 34.0
        v_['A39'] = 36.0
        v_['A40'] = 36.0
        v_['A41'] = 38.0
        v_['A42'] = 38.0
        v_['A43'] = 40.0
        v_['A44'] = 42.0
        v_['B1'] = 0.49
        v_['B2'] = 0.49
        v_['B3'] = 0.48
        v_['B4'] = 0.47
        v_['B5'] = 0.48
        v_['B6'] = 0.47
        v_['B7'] = 0.46
        v_['B8'] = 0.46
        v_['B9'] = 0.45
        v_['B10'] = 0.43
        v_['B11'] = 0.45
        v_['B12'] = 0.43
        v_['B13'] = 0.43
        v_['B14'] = 0.44
        v_['B15'] = 0.43
        v_['B16'] = 0.43
        v_['B17'] = 0.46
        v_['B18'] = 0.45
        v_['B19'] = 0.42
        v_['B20'] = 0.42
        v_['B21'] = 0.43
        v_['B22'] = 0.41
        v_['B23'] = 0.41
        v_['B24'] = 0.40
        v_['B25'] = 0.42
        v_['B26'] = 0.40
        v_['B27'] = 0.40
        v_['B28'] = 0.41
        v_['B29'] = 0.40
        v_['B30'] = 0.41
        v_['B31'] = 0.41
        v_['B32'] = 0.40
        v_['B33'] = 0.40
        v_['B34'] = 0.40
        v_['B35'] = 0.38
        v_['B36'] = 0.41
        v_['B37'] = 0.40
        v_['B38'] = 0.40
        v_['B39'] = 0.41
        v_['B40'] = 0.38
        v_['B41'] = 0.40
        v_['B42'] = 0.40
        v_['B43'] = 0.39
        v_['B44'] = 0.39
        v_['1'] = 1
        v_['44'] = 44
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(0.49))
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
        self.gconst = arrset(self.gconst,ig_['CON1'],float(0.09))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.4
        self.xlower[ix_['X2']] = -4.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(0.42)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(0.42)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(5.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(5.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'AA')
        elftp = loaset(elftp,it,1,'BB')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['44'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eOBSQ')
            ielftype = arrset(ielftype,ie,iet_["eOBSQ"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='AA')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(I)]))
            posep = np.where(elftp[ielftype[ie]]=='BB')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['B'+str(I)]))
        ename = 'PR'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['44'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CON1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PR'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                 0.02845966
# LO SOLTN                 0.03063791
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSQR2-AN-2-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eOBSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        AM8 = self.elpar[iel_][0]-8.0
        CMV1 = 0.49-EV_[0,0]
        E = np.exp(-EV_[1,0]*AM8)
        DED2 = -AM8*E
        R = self.elpar[iel_][1]-EV_[0,0]-CMV1*E
        DRD1 = E-1.0
        DRD2 = -CMV1*DED2
        D2RD22 = -CMV1*AM8*AM8*E
        f_   = R*R
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*R*DRD1
            g_[1] = 2.0*R*DRD2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*DRD1*DRD1
                H_[0,1] = 2.0*(DRD2*DRD1+R*DED2)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*(DRD2*DRD2+R*D2RD22)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

