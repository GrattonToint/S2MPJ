from s2mpjlib import *
class  HS117(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS117
#    *********
# 
#    Source: problem 117 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-COQR2-AN-15-5"
# 
#    Number of constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS117'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 5
        v_['M'] = 10
        v_['M+N'] = v_['M']+v_['N']
        v_['1'] = 1
        v_['E1'] = -15.0
        v_['E2'] = -27.0
        v_['E3'] = -36.0
        v_['E4'] = -18.0
        v_['E5'] = -12.0
        v_['C1,1'] = 30.0
        v_['C2,1'] = -20.0
        v_['C3,1'] = -10.0
        v_['C4,1'] = 32.0
        v_['C5,1'] = -10.0
        v_['C1,2'] = -20.0
        v_['C2,2'] = 39.0
        v_['C3,2'] = -6.0
        v_['C4,2'] = -31.0
        v_['C5,2'] = 32.0
        v_['C1,3'] = -10.0
        v_['C2,3'] = -6.0
        v_['C3,3'] = 10.0
        v_['C4,3'] = -6.0
        v_['C5,3'] = -10.0
        v_['C1,4'] = 32.0
        v_['C2,4'] = -31.0
        v_['C3,4'] = -6.0
        v_['C4,4'] = 39.0
        v_['C5,4'] = -20.0
        v_['C1,5'] = -10.0
        v_['C2,5'] = 32.0
        v_['C3,5'] = -10.0
        v_['C4,5'] = -20.0
        v_['C5,5'] = 30.0
        v_['D1'] = 4.0
        v_['D2'] = 8.0
        v_['D3'] = 10.0
        v_['D4'] = 6.0
        v_['D5'] = 2.0
        v_['A1,1'] = -16.0
        v_['A2,1'] = 0.0
        v_['A3,1'] = -3.5
        v_['A4,1'] = 0.0
        v_['A5,1'] = 0.0
        v_['A6,1'] = 2.0
        v_['A7,1'] = -1.0
        v_['A8,1'] = -1.0
        v_['A9,1'] = 1.0
        v_['A10,1'] = 1.0
        v_['A1,2'] = 2.0
        v_['A2,2'] = -2.0
        v_['A3,2'] = 0.0
        v_['A4,2'] = -2.0
        v_['A5,2'] = -9.0
        v_['A6,2'] = 0.0
        v_['A7,2'] = -1.0
        v_['A8,2'] = -2.0
        v_['A9,2'] = 2.0
        v_['A10,2'] = 1.0
        v_['A1,3'] = 0.0
        v_['A2,3'] = 0.0
        v_['A3,3'] = 2.0
        v_['A4,3'] = 0.0
        v_['A5,3'] = -2.0
        v_['A6,3'] = -4.0
        v_['A7,3'] = -1.0
        v_['A8,3'] = -3.0
        v_['A9,3'] = 3.0
        v_['A10,3'] = 1.0
        v_['A1,4'] = 1.0
        v_['A2,4'] = 4.0
        v_['A3,4'] = 0.0
        v_['A4,4'] = -4.0
        v_['A5,4'] = 1.0
        v_['A6,4'] = 0.0
        v_['A7,4'] = -1.0
        v_['A8,4'] = -2.0
        v_['A9,4'] = 4.0
        v_['A10,4'] = 1.0
        v_['A1,5'] = 0.0
        v_['A2,5'] = 2.0
        v_['A3,5'] = 0.0
        v_['A4,5'] = -1.0
        v_['A5,5'] = -2.8
        v_['A6,5'] = 0.0
        v_['A7,5'] = -1.0
        v_['A8,5'] = -1.0
        v_['A9,5'] = 5.0
        v_['A10,5'] = 1.0
        v_['B1'] = -40.0
        v_['B2'] = -2.0
        v_['B3'] = -0.25
        v_['B4'] = -4.0
        v_['B5'] = -4.0
        v_['B6'] = -1.0
        v_['B7'] = -40.0
        v_['B8'] = -60.0
        v_['B9'] = 5.0
        v_['B10'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['M+N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            v_['-BJ'] = -1.0*v_['B'+str(J)]
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(J)]])
            valA = np.append(valA,float(v_['-BJ']))
        for J in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['N'])+1):
                v_['M+K'] = v_['M']+K
                v_['2CKJ'] = 2.0*v_['C'+str(K)+','+str(J)]
                [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['M+K']))]])
                valA = np.append(valA,float(v_['2CKJ']))
            for K in range(int(v_['1']),int(v_['M'])+1):
                v_['-AKJ'] = -1.0*v_['A'+str(K)+','+str(J)]
                [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(K)]])
                valA = np.append(valA,float(v_['-AKJ']))
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
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['-EJ'] = -1.0*v_['E'+str(J)]
            self.gconst = arrset(self.gconst,ig_['C'+str(J)],float(v_['-EJ']))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.001))
        self.y0 = np.full((self.m,1),float(0.001))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(60.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X7']),float(60.0)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQUARE', iet_)
        elftv = loaset(elftv,it,0,'XJ')
        [it,iet_,_] = s2mpj_ii( 'eCUBE', iet_)
        elftv = loaset(elftv,it,0,'XJ')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'XI')
        elftv = loaset(elftv,it,1,'XJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['M+J'] = v_['M']+J
            ename = '2D'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCUBE')
            ielftype = arrset(ielftype,ie,iet_["eCUBE"])
            vname = 'X'+str(int(v_['M+J']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.001))
            posev = np.where(elftv[ielftype[ie]]=='XJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = '3D'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQUARE')
            ielftype = arrset(ielftype,ie,iet_["eSQUARE"])
            vname = 'X'+str(int(v_['M+J']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.001))
            posev = np.where(elftv[ielftype[ie]]=='XJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for K in range(int(v_['1']),int(v_['N'])+1):
                v_['M+K'] = v_['M']+K
                ename = 'C'+str(K)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'ePROD')
                ielftype = arrset(ielftype,ie,iet_["ePROD"])
                vname = 'X'+str(int(v_['M+K']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.001))
                posev = np.where(elftv[ielftype[ie]]=='XI')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['M+J']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.001))
                posev = np.where(elftv[ielftype[ie]]=='XJ')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['2DJ'] = 2.0*v_['D'+str(J)]
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['2D'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['2DJ']))
            v_['3DJ'] = 3.0*v_['D'+str(J)]
            ig = ig_['C'+str(J)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['3D'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['3DJ']))
            for K in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['OBJ']
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(K)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['C'+str(K)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               32.34867897
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COQR2-AN-15-5"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQUARE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e+0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCUBE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0e+0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0e+0
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

