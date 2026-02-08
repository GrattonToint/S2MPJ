from s2mpjlib import *
class  HS85(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS85
#    *********
# 
#    The problem is to optimize the net profit of an hypothetical
#    wood-pulp plant. The constraints include the usual material
#    and energy balances as well as several empirical equations.
# 
#    Source: problem 85 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#               Python coding: Cunxin Huang, 2025.
# 
#    classification = "C-COOI2-MN-5-21"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS85'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)
        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N']   = 5
        v_['A2']  = 17.505
        v_['A3']  = 11.275
        v_['A4']  = 214.228
        v_['A5']  = 7.458
        v_['A6']  = 0.961
        v_['A7']  = 1.612
        v_['A8']  = 0.146
        v_['A9']  = 107.99
        v_['A10'] = 922.693
        v_['A11'] = 926.832
        v_['A12'] = 18.766
        v_['A13'] = 1072.163
        v_['A14'] = 8961.448
        v_['A15'] = 0.063
        v_['A16'] = 71084.33
        v_['A17'] = 2802713.0
        v_['B2']  = 1053.6667
        v_['B3']  = 35.03
        v_['B4']  = 665.585
        v_['B5']  = 584.463
        v_['B6']  = 265.916
        v_['B7']  = 7.046
        v_['B8']  = 0.222
        v_['B9']  = 273.366
        v_['B10'] = 1286.105
        v_['B11'] = 1444.046
        v_['B12'] = 537.141
        v_['B13'] = 3247.039
        v_['B14'] = 26844.086
        v_['B15'] = 0.386
        v_['B16'] = 140000.0
        v_['B17'] = 12146108.0
        v_['1']   = 1
        v_['2']   = 2
        v_['17']  = 17
        v_['19']  = 19
        v_['20']  = 20
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars     = np.array([])
        binvars     = np.array([])
        irA         = np.array([],dtype=int)
        icA         = np.array([],dtype=int)
        valA        = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CON0',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON0')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(1.5))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['20'])+1):
            [ig,ig_,_] = s2mpj_ii('CON'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CON'+str(I))
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
        self.gconst = arrset(self.gconst,ig_['OBJ'],float(0.1365))
        self.gconst = arrset(self.gconst,ig_['CON1'],float(213.1))
        for I in range(int(v_['2']),int(v_['17'])+1):
            self.gconst = arrset(self.gconst,ig_['CON'+str(I)],float(v_['A'+str(I)]))
        self.gconst = arrset(self.gconst,ig_['CON19'],float(-21.0))
        self.gconst = arrset(self.gconst,ig_['CON20'],float(110.6))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        for I in range(int(v_['2']),int(v_['17'])+1):
            v_['DIF'] = v_['B'+str(I)]-v_['A'+str(I)]
            grange = arrset(grange,ig_['CON'+str(I)],float(v_['DIF']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 7.044148e+2
        self.xupper[ix_['X1']] = 9.063855e+2
        self.xlower[ix_['X2']] = 6.86e+1
        self.xupper[ix_['X2']] = 2.8888e+2
        self.xupper[ix_['X3']] = 1.3475e+2
        self.xlower[ix_['X4']] = 1.930e+2
        self.xupper[ix_['X4']] = 2.870966e+2
        self.xlower[ix_['X5']] = 2.50e+1
        self.xupper[ix_['X5']] = 8.41988e+1
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(9.0e+2)
        else:
            self.y0  = arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(9.0e+2))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(8.0e+1)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(8.0e+1))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(1.15e+2)
        else:
            self.y0  = arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(1.15e+2))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(2.67e+2)
        else:
            self.y0  = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(2.67e+2))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(2.7e+1)
        else:
            self.y0  = arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(2.7e+1))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eY',iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftp = []
        elftp = loaset(elftp,it,0,'PI')
        [it,iet_,_] = s2mpj_ii( 'eC',iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftp = loaset(elftp,it,0,'PI')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['19'])+1):
            v_['PI'] = float(I)
            v_['PI'] = 0.01+v_['PI']
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eC')
            ielftype = arrset(ielftype,ie,iet_["eC"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='PI')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PI']))
        for I in range(int(v_['1']),int(v_['20'])+1):
            v_['PI'] = float(I)
            v_['PI'] = 0.01+v_['PI']
            ename = 'Y'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eY')
            ielftype = arrset(ielftype,ie,iet_["eY"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='PI')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PI']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc          = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y17'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-5.843e-7))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y14'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.17e-4))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y13'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.358e-5))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y16'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.502e-6))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.0321))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['Y5'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.00423))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0e-4))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C19'])
        self.grelw = loaset(self.grelw,ig,posel,float(37.48))
        for I in range(int(v_['1']),int(v_['20'])+1):
            ig = ig_['CON'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['Y'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1.90513375
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nge)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0]
        self.pbclass   = "C-COOI2-MN-5-21"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

#    @staticmethod
#    def e_globs(self):

#        import numpy as np
#        self.efpar = np.array([]);
#        self.efpar = arrset( self.efpar,0,0)
#        return pbm
        
    @staticmethod
    def eY(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        I    = int(self.elpar[iel_][0]) - 1
        C, Y = self.extfunc(self,EV_)
        f_   = Y[0,I]
        if nargout>1:
            dim   = len(EV_)
            g_    = np.zeros(dim)
            g_[0] = Y[1,I]
            g_[1] = Y[2,I]
            g_[2] = Y[3,I]
            g_[3] = Y[4,I]
            g_[4] = Y[5,I]
            if nargout>2:
                H_      = np.zeros((5,5))
                H_[0,0] = Y[6,I]
                H_[0,1] = Y[7,I]
                H_[1,0] = H_[0,1]
                H_[1,1] = Y[8,I]
                H_[0,2] = Y[9,I]
                H_[2,0] = H_[0,2]
                H_[1,2] = Y[10,I]
                H_[2,1] = H_[1,2]
                H_[2,2] = Y[11,I]
                H_[0,3] = Y[12,I]
                H_[3,0] = H_[0,3]
                H_[1,3] = Y[13,I]
                H_[3,1] = H_[1,3]
                H_[2,3] = Y[14,I]
                H_[3,2] = H_[2,3]
                H_[3,3] = Y[15,I]
                H_[0,4] = Y[16,I]
                H_[4,0] = H_[0,4]
                H_[1,4] = Y[17,I]
                H_[4,1] = H_[1,4]
                H_[2,4] = Y[18,I]
                H_[4,2] = H_[2,4]
                H_[3,4] = Y[19,I]
                H_[4,3] = H_[3,4]
                H_[4,4] = Y[20,I]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        I    = int(self.elpar[iel_][0]) - 1
        C, Y = self.extfunc(self,EV_)
        f_   = C[0,I]
        if nargout>1:
            dim   = len(EV_)
            g_    = np.zeros(dim)
            g_[0] = C[1,I]
            g_[1] = C[2,I]
            g_[2] = C[3,I]
            g_[3] = C[4,I]
            g_[4] = C[5,I]
            if nargout>2:
                H_      = np.zeros((5,5))
                H_[0,0] = C[6,I]
                H_[0,1] = C[7,I]
                H_[1,0] = H_[0,1]
                H_[1,1] = C[8,I]
                H_[0,2] = C[9,I]
                H_[2,0] = H_[0,2]
                H_[1,2] = C[10,I]
                H_[2,1] = H_[1,2]
                H_[2,2] = C[11,I]
                H_[0,3] = C[12,I]
                H_[3,0] = H_[0,3]
                H_[1,3] = C[13,I]
                H_[3,1] = H_[1,3]
                H_[2,3] = C[14,I]
                H_[3,2] = H_[2,3]
                H_[3,3] = C[15,I]
                H_[0,4] = C[16,I]
                H_[4,0] = H_[0,4]
                H_[1,4] = C[17,I]
                H_[4,1] = H_[1,4]
                H_[2,4] = C[18,I]
                H_[4,2] = H_[2,4]
                H_[3,4] = C[19,I]
                H_[4,3] = H_[3,4]
                H_[4,4] = C[20,I]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_
        
    @staticmethod
    def extfunc(self,X):
        # A translation of the Fortran code present in the SIF file.
        import numpy as np

        X1 = X[0].item()
        X2 = X[1].item()
        X3 = X[2].item()
        X4 = X[3].item()
        X5 = X[4].item()

        C = np.zeros((21,19))
        Y = np.zeros((21,20))
        
        # Function Y1.

        A      = 4.16e+1
        Y[0,0] = X2 + X3 + A
        Y[2,0] = 1.0
        Y[3,0] = 1.0

        # Function C1.

        A      = 2.4e-2
        B      = 4.62
        C[0,0] = A*X4 - B
        C[4,0] = A

        # Function Y2.

        A       = 1.2e+1
        B       = 1.25e+1
        Y[0,1]  = A + B/C[0,0]
        Y[4,1]  = - B*C[4,0]/C[0,0]**2
        Y[15,1] = - B*C[15,0]/C[0,0]**2 + 2.0*B*C[4,0]**2/C[0,0]**3

        # Function C2.

        A       = 3.535e-4
        B       = 5.311e-1
        D       = 8.705e-2
        C[0,1]  = A*X1**2 + B*X1 + D*X1*Y[0,1]
        C[1,1]  = 2.0*A*X1 + B + D*Y[0,1]
        C[4,1]  = D*X1*Y[4,1]
        C[6,1]  = 2.0*A
        C[12,1] = D*Y[4,1]
        C[15,1] = D*X1*Y[15,1]

        # Function C3.

        A       = 5.2e-2
        B       = 7.8e+1
        D       = 2.377e-3
        C[0,2]  = A*X1 + B + D*X1*Y[0,1]
        C[1,2]  = A + D*Y[0,1]
        C[4,2]  = D*X1*Y[4,1]
        C[12,2] = D*Y[4,1]
        C[15,2] = D*X1*Y[15,1]

        # Function Y3.

        Y[0,2]  = C[0,1]/C[0,2]
        Y[1,2]  = C[1,1]/C[0,2] - C[0,1]*C[1,2]/C[0,2]**2
        Y[4,2]  = C[4,1]/C[0,2] - C[0,1]*C[4,2]/C[0,2]**2
        Y[6,2]  = ( C[6,1]/C[0,2] - 2.0*C[1,1]*C[1,2]/C[0,2]**2  
                  - C[0,1]*C[6,2]/C[0,2]**2 + 2.0*C[0,1]*C[1,2]**2/C[0,2]**3 )
        Y[12,2] = ( C[12,1]/C[0,2] -  C[1,1]*C[4,2]/C[0,2]**2 - C[4,1]*C[1,2]/C[0,2]**2
                  - C[0,1]*C[12,2]/C[0,2]**2 + 2.0*C[0,1]*C[1,2]*C[4,2]/C[0,2]**3 )
        Y[15,2] = ( C[15,1]/C[0,2] - 2.0*C[4,1]*C[4,2]/C[0,2]**2 
                  - C[0,1]*C[15,2]/C[0,2]**2 + 2.0*C[0,1]*C[4,2]**2/C[0,2]**3 )

        # Function Y4.

        A       = 1.9e+1
        Y[0,3]  = A*Y[0,2]
        Y[1,3]  = A*Y[1,2]
        Y[4,3]  = A*Y[4,2]
        Y[6,3]  = A*Y[6,2]
        Y[12,3] = A*Y[12,2]
        Y[15,3] = A*Y[15,2]

        # Function C4.

        A       = 4.782e-2
        B       = 1.956e-1
        D       = 6.376e-1
        E       = 1.594
        C[0,3]  = A*( X1 - Y[0,2]) +  B*( ( X1 - Y[0,2])**2)/X2 + D*Y[0,3] + E*Y[0,2]
        C[1,3]  = A*( 1.0 - Y[1,2]) + 2.0*B*( X1 - Y[0,2])*( 1.0 - Y[1,2])/X2 + D*Y[1,3] + E*Y[1,2]
        C[2,3]  = - B*( ( X1 - Y[0,2])**2)/X2**2
        C[4,3]  = A*( - Y[4,2]) + 2.0*B*( X1 - Y[0,2])*( - Y[4,2])/X2 + D*Y[4,3] + E*Y[4,2]
        C[6,3]  = ( A*( - Y[6,2]) + 2.0*B*( 1.0 - Y[1,2])**2/X2 + 2.0*B*( X1 - Y[0,2])*( - Y[6,2])/X2
                  + D*Y[6,3] + E*Y[6,2] )
        C[7,3]  = - 2.0*B*( X1 - Y[0,2])*( 1.0 - Y[1,2])/X2**2
        C[8,3]  = 2.0*B*( ( X1 - Y[0,2])**2)/X2**3
        C[12,3] = ( A*( - Y[12,2]) + 2.0*B*( - Y[4,2])*( 1.0 - Y[1,2])/X2
                  + 2.0*B*( X1 - Y[0,2])*( - Y[12,2])/X2 + D*Y[12,3] + E*Y[12,2] )
        C[13,3] = - 2.0*B*( X1 - Y[0,2])*( - Y[4,2])/X2**2
        C[15,3] = ( A*( - Y[15,2]) + 2.0*B*( - Y[4,2])**2/X2 + 2.0*B*( X1 - Y[0,2])*( - Y[15,2])/X2
                  + D*Y[15,3] + E*Y[15,2] )
        # Function C5.

        A       = 1.0e+2
        C[0,4]  = A*X2
        C[2,4]  = A

        # Function C6.

        C[0,5]  = X1  - Y[0,2] - Y[0,3]
        C[1,5]  = 1.0 - Y[1,2] - Y[1,3]
        C[4,5]  = - Y[4,2] - Y[4,3]
        C[6,5]  = - Y[6,2] - Y[6,3]
        C[12,5] = - Y[12,2] - Y[12,3]
        C[15,5] = - Y[15,2] - Y[15,3]

        # Function C7.

        A       = 9.5e-1
        C[0,6]  = A - C[0,3]/C[0,4]
        C[1,6]  = - C[1,3]/C[0,4]
        C[2,6]  = - C[2,3]/C[0,4] + C[0,3]*C[2,4]/C[0,4]**2
        C[4,6]  = - C[4,3]/C[0,4]
        C[6,6]  = - C[6,3]/C[0,4]
        C[7,6]  = - C[7,3]/C[0,4] + C[1,3]*C[2,4]/C[0,4]**2
        C[8,6]  = ( - C[8,3]/C[0,4] + 2.0*C[2,3]*C[2,4]/C[0,4]**2 + C[0,3]*C[8,4]/C[0,4]**2
                  - 2.0*C[0,3]*C[2,4]**2/C[0,4]**3 )
        C[12,6] = - C[12,3]/C[0,4]
        C[13,6] = - C[13,3]/C[0,4] + C[4,3]*C[2,4]/C[0,4]**2 + C[0,3]*C[13,4]/C[0,4]**2
        C[15,6] = - C[15,3]/C[0,4]

        # Function Y5.

        Y[0,4]  = C[0,5]*C[0,6]
        Y[1,4]  = C[1,5]*C[0,6] + C[0,5]*C[1,6]
        Y[2,4]  = C[0,5]*C[2,6]
        Y[4,4]  = C[4,5]*C[0,6] + C[0,5]*C[4,6]
        Y[6,4]  = C[6,5]*C[0,6] + 2.0*C[1,5]*C[1,6] + C[0,5]*C[6,6]
        Y[7,4]  = C[1,5]*C[2,6] + C[0,5]*C[7,6]
        Y[8,4]  = C[0,5]*C[8,6]
        Y[12,4] = C[12,5]*C[0,6] + C[4,5]*C[1,6] + C[1,5]*C[4,6] + C[0,5]*C[12,6]
        Y[13,4] = C[4,5] *C[2,6] + C[0,5]*C[13,6]
        Y[15,4] = C[15,5]*C[0,6] + 2.0*C[4,5]*C[4,6] + C[0,5]*C[15,6]

        # Function Y6.

        Y[0,5]  = X1 - Y[0,2] - Y[0,3] - Y[0,4]
        Y[1,5]  = 1.0 - Y[1,2] - Y[1,3] - Y[1,4]
        Y[2,5]  = - Y[2,4]
        Y[4,5]  = - Y[4,2] - Y[4,3] - Y[4,4]
        Y[6,5]  = - Y[6,2] - Y[6,3] - Y[6,4]
        Y[7,5]  = - Y[7,4]
        Y[8,5]  = - Y[8,4]
        Y[12,5] = - Y[12,2] - Y[12,3] - Y[12,4]
        Y[13,5] = - Y[13,4]
        Y[15,5] = - Y[15,2] - Y[15,3] - Y[15,4]

        # Function C8.

        A       = 9.95e-1
        C[0,7]  = A*( Y[0,3] + Y[0,4])
        C[1,7]  = A*( Y[1,3] + Y[1,4])
        C[2,7]  = A*Y[2,4]
        C[4,7]  = A*( Y[4,3] + Y[4,4])
        C[6,7]  = A*( Y[6,3] + Y[6,4])
        C[7,7]  = A*Y[7,4]
        C[8,7]  = A*Y[8,4]
        C[12,7] = A*( Y[12,3] + Y[12,4])
        C[13,7] = A*Y[13,4]
        C[15,7] = A*( Y[15,3] + Y[15,4])

        # Function Y7.

        Y[0,6]  = C[0,7]/Y[0,0]
        Y[1,6]  = C[1,7]/Y[0,0]
        Y[2,6]  = C[2,7]/Y[0,0] - C[0,7]*Y[2,0]/Y[0,0]**2
        Y[3,6]  = - C[0,7]*Y[3,0]/Y[0,0]**2
        Y[4,6]  = C[4,7]/Y[0,0]
        Y[6,6]  = C[6,7]/Y[0,0]
        Y[7,6]  = C[7,7]/Y[0,0] - C[1,7]*Y[2,0]/Y[0,0]**2
        Y[8,6]  = C[8,7]/Y[0,0] - 2.0*C[2,7]*Y[2,0]/Y[0,0]**2 + 2.0*C[0,7]*Y[2,0]**2/Y[0,0]**3
        Y[9,6]  = - C[1,7]*Y[3,0]/Y[0,0]**2
        Y[10,6] = ( C[10,7]/Y[0,0] -  C[2,7]*Y[3,0]/Y[0,0]**2 - C[3,7]*Y[2,0]/Y[0,0]**2
                  + 2.0*C[0,7]*Y[2,0]*Y[3,0]/Y[0,0]**3 )
        Y[11,6] = - C[3,7]*Y[3,0]/Y[0,0]**2 + 2.0*C[0,7]*Y[3,0]**2/Y[0,0]**3
        Y[12,6] = C[12,7]/Y[0,0]
        Y[13,6] = C[13,7]/Y[0,0] - C[4,7]*Y[2,0]/Y[0,0]**2
        Y[14,6] = - C[4,7]*Y[3,0]/Y[0,0]**2 + 2.0*C[0,7]*Y[3,0]*Y[4,0]/Y[0,0]**3
        Y[15,6] = C[15,7]/Y[0,0]

        # Function Y8.

        A       = 3.798e+3
        Y[0,7]  = C[0,7] /A
        Y[1,7]  = C[1,7] /A
        Y[2,7]  = C[2,7] /A
        Y[4,7]  = C[4,7] /A
        Y[6,7]  = C[6,7] /A
        Y[7,7]  = C[7,7] /A
        Y[8,7]  = C[8,7] /A
        Y[12,7] = C[12,7]/A
        Y[13,7] = C[13,7]/A
        Y[15,7] = C[15,7]/A

        # Function C9.

        A       = 6.63e-2
        B       = 3.153e-1
        C[0,8]  = Y[0,6] - A*Y[0,6]/Y[0,7] - B
        C[1,8]  = Y[1,6] - A*Y[1,6]/Y[0,7] + A*Y[0,6]*Y[1,7]/Y[0,7]**2
        C[2,8]  = Y[2,6] - A*Y[2,6]/Y[0,7] + A*Y[0,6]*Y[2,7]/Y[0,7]**2
        C[3,8]  = Y[3,6] - A*Y[3,6]/Y[0,7] + A*Y[0,6]*Y[3,7]/Y[0,7]**2
        C[4,8]  = Y[4,6] - A*Y[4,6]/Y[0,7] + A*Y[0,6]*Y[4,7]/Y[0,7]**2
        C[6,8]  = ( Y[6,6] - A*Y[6,6]/Y[0,7] + 2.0*A*Y[1,6]*Y[1,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[6,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[1,7]**2/Y[0,7]**3 )
        C[7,8]  = ( Y[7,6] - A*Y[7,6]/Y[0,7] + A*Y[1,6]*Y[2,7]/ Y[0,7]**2
                  + A*Y[2,6]*Y[1,7]/Y[0,7]**2 + A*Y[0,6]*Y[7,7]/Y[0,7]**2
                  - 2.0*A*Y[0,6]*Y[1,7]*Y[2,7]/Y[0,7]**3 )
        C[8,8]  = ( Y[8,6] - A*Y[8,6]/Y[0,7] +  2.0*A*Y[2,6]*Y[2,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[8,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[2,7]**2/Y[0,7]**3 )
        C[9,8]  = ( Y[9,6] - A*Y[9,6]/Y[0,7]
                  + A*Y[1,6]*Y[3,7]/ Y[0,7]**2 +  A*Y[3,6]*Y[1,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[9,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[1,7]*Y[3,7]/Y[0,7]**3 )
        C[10,8] = ( Y[10,6] - A*Y[10,6]/Y[0,7]
                  + A*Y[2,6]*Y[3,7]/ Y[0,7]**2 + A*Y[3,6]*Y[2,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[10,7]/Y[0,7]**2 -2.0*A*Y[0,6]*Y[2,7]*Y[3,7]/Y[0,7]**3 )
        C[11,8] = ( Y[11,6] - A*Y[11,6]/Y[0,7] + 2.0*A*Y[3,6]*Y[3,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[11,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[3,7]**2/Y[0,7]**3 )
        C[12,8] = ( Y[12,6] - A*Y[12,6]/Y[0,7]
                  + A*Y[1,6]*Y[4,7]/ Y[0,7]**2 + A*Y[4,6]*Y[1,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[12,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[1,7]*Y[4,7]/Y[0,7]**3 )
        C[13,8] = ( Y[13,6] - A*Y[13,6]/Y[0,7]
                  + A*Y[2,6]*Y[4,7]/ Y[0,7]**2 + A*Y[4,6]*Y[2,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[13,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[2,7]*Y[4,7]/Y[0,7]**3 )
        C[14,8] = ( Y[14,6] - A*Y[14,6]/Y[0,7]
                  + A*Y[3,6]*Y[4,7]/ Y[0,7]**2 + A*Y[4,6]*Y[3,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[14,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[3,7]*Y[4,7]/Y[0,7]**3 )
        C[15,8] = ( Y[15,6] - A*Y[15,6]/Y[0,7] + 2.0*A*Y[4,6]*Y[4,7]/Y[0,7]**2
                  + A*Y[0,6]*Y[15,7]/Y[0,7]**2 - 2.0*A*Y[0,6]*Y[4,7]**2/Y[0,7]**3 )

        # Function Y9.

        A       = 9.682e+1
        B       = 3.21e-1
        Y[0,8]  = A/C[0,8] + B*Y[0,0]
        Y[1,8]  = - A*C[1,8]/C[0,8]**2
        Y[2,8]  = - A*C[2,8]/C[0,8]**2 + B*Y[2,0]
        Y[3,8]  = - A*C[3,8]/C[0,8]**2 + B*Y[3,0]
        Y[4,8]  = - A*C[4,8]/C[0,8]**2
        Y[6,8]  = - A*C[6,8]/C[0,8]**2  + 2.0*A*C[1,8]**2/C[0,8]**3
        Y[7,8]  = - A*C[7,8]/C[0,8]**2  + 2.0*A*C[1,8]*C[2,8]/C[0,8]**3
        Y[8,8]  = - A*C[8,8]/C[0,8]**2  + 2.0*A*C[2,8]**2/C[0,8]**3
        Y[9,8]  = - A*C[9,8]/C[0,8]**2 + 2.0*A*C[1,8]*C[3,8]/C[0,8]**3
        Y[10,8] = - A*C[10,8]/C[0,8]**2 + 2.0*A*C[2,8]*C[3,8]/C[0,8]**3
        Y[11,8] = - A*C[11,8]/C[0,8]**2 + 2.0*A*C[3,8]**2/C[0,8]**3
        Y[12,8] = - A*C[12,8]/C[0,8]**2 + 2.0*A*C[1,8]*C[4,8]/C[0,8]**3
        Y[13,8] = - A*C[13,8]/C[0,8]**2 + 2.0*A*C[2,8]*C[4,8]/C[0,8]**3
        Y[14,8] = - A*C[14,8]/C[0,8]**2 + 2.0*A*C[3,8]*C[4,8]/C[0,8]**3
        Y[15,8] = - A*C[15,8]/C[0,8]**2 + 2.0*A*C[4,8]**2/C[0,8]**3

        # Function Y10.

        A       = 2.29
        B       = 1.258
        D       = 1.29
        E       = 1.71
        Y[0,9]  = A*Y[0,2]  + B*Y[0,3]  + D*Y[0,4] + E*Y[0,5]
        Y[1,9]  = A*Y[1,2]  + B*Y[1,3]  + D*Y[1,4] + E*Y[1,5]
        Y[2,9]  = D*Y[2,4]  + E*Y[2,5]
        Y[4,9]  = A*Y[4,2]  + B*Y[4,3]  + D*Y[4,4] + E*Y[4,5]
        Y[6,9]  = A*Y[6,2]  + B*Y[6,3]  + D*Y[6,4] + E*Y[6,5]
        Y[7,9]  = D*Y[7,4]  + E*Y[7,5]
        Y[8,9]  = D*Y[8,4]  + E*Y[8,5]
        Y[12,9] = A*Y[12,2] + B*Y[12,3] + D*Y[12,4] + E*Y[12,5]
        Y[13,9] = D*Y[13,4] + E*Y[13,5]
        Y[15,9] = A*Y[15,2] + B*Y[15,3] + D*Y[15,4] + E*Y[15,5]

        # Function Y11.

        A        = 1.71
        B        = 5.8e-1
        D        = 4.52e-1
        Y[0,10]  = A*X1 + B*Y[0,2] - D*Y[0,3]
        Y[1,10]  = A + B*Y[1,2] - D*Y[1,3]
        Y[4,10]  = B*Y[4,2]  - D*Y[4,3]
        Y[6,10]  = B*Y[6,2]  - D*Y[6,3]
        Y[12,10] = B*Y[12,2] - D*Y[12,3]
        Y[15,10] = B*Y[15,2] - D*Y[15,3]

        # Function C10.

        A       = 1.23e+1
        B       = 7.523e+2
        C[0,9]  = A/B

        # Function C11.

        A        = 1.74125
        C[0,10]  = A*X1*Y[0,1]
        C[1,10]  = A*Y[0,1]
        C[4,10]  = A*X1*Y[4,1]
        C[12,10] = A*Y[4,1]
        C[15,10] = A*X1*Y[15,1]

        # Function C12.

        A        = 9.995e-1
        B        = 1.998e+3
        C[0,11]  = A*Y[0,9] + B
        C[1,11]  = A*Y[1,9]
        C[2,11]  = A*Y[2,9]
        C[4,11]  = A*Y[4,9]
        C[6,11]  = A*Y[6,9]
        C[7,11]  = A*Y[7,9]
        C[8,11]  = A*Y[8,9]
        C[12,11] = A*Y[12,9]
        C[13,11] = A*Y[13,9]
        C[15,11] = A*Y[15,9]

        # Function Y12.

        Y[0,11]  = C[0,9]*X1 + C[0,10]/C[0,11]
        Y[1,11]  = C[0,9] + C[1,10]/C[0,11] - C[0,10]*C[1,11]/C[0,11]**2
        Y[2,11]  = - C[0,10]*C[2,11]/C[0,11]**2
        Y[4,11]  = C[4,10]/C[0,11] - C[0,10]*C[4,11]/C[0,11]**2
        Y[6,11]  = ( C[6,10]/C[0,11] - 2.0*C[1,10]*C[1,11]/C[0,11]**2
                   - C[0,10]*C[6,11]/C[0,11]**2 + 2.0*C[0,10]*C[1,11]**2/C[0,11]**3 )
        Y[7,11]  = ( C[1,10]*C[2,11]/C[0,11]**2 -C[0,10]*C[7,11]/C[0,11]**2
                   + 2.0*C[0,10]*C[1,11]*C[2,11]/C[0,11]**3 )
        Y[8,11]  = - C[0,10]*C[8,11]/C[0,11]**2 + 2.0*C[0,10]*C[2,11]**2/C[0,11]**3
        Y[12,11] = ( C[12,10]/C[0,11] - C[1,10]*C[4,11]/C[0,11]**2
                   - C[4,10]*C[1,11]/C[0,11]**2 - C[0,10]*C[12,11]/C[0,11]**2
                   + 2.0*C[0,10]*C[1,11]*C[4,11]/C[0,11]**3 )
        Y[13,11] = (- C[4,10]*C[2,11]/C[0,11]**2 -  C[0,10]*C[13,11]/C[0,11]**2
                   + 2.0*C[0,10]*C[2,11]*C[4,11]/C[0,11]**3 )
        Y[15,11] = ( C[15,10]/C[0,11] - 2.0*C[4,10]*C[4,11]/C[0,11]**2
                   - C[0,10]*C[15,11]/C[0,11]**2 + 2.0*C[0,10]*C[4,11]**2/C[0,11]**3 )

        # Function Y13.

        A        = 1.75
        Y[0,12]  = C[0,11] - A*Y[0,1]
        Y[1,12]  = C[1,11]
        Y[2,12]  = C[2,11]
        Y[4,12]  = C[4,11] - A*Y[4,1]
        Y[6,12]  = C[6,11]
        Y[7,12]  = C[7,11]
        Y[8,12]  = C[8,11]
        Y[12,12] = C[12,11]
        Y[13,12] = C[13,11]
        Y[15,12] = C[15,11] - A*Y[15,1]

        # Function Y14.

        A        = 3.623e+3
        B        = 6.44e+1
        D        = 1.46312e+5
        F        = 5.84e+1
        Y[0,13]  = A + B*X2 + F*X3 + D/( Y[0,8] + X5 )
        Y[1,13]  = - D*Y[1,8]/( Y[0,8] + X5 )**2
        Y[2,13]  = B - D*Y[2,8]/( Y[0,8] + X5 )**2
        Y[3,13]  = F - D*Y[3,8]/( Y[0,8] + X5 )**2
        Y[4,13]  = - D*Y[4,8]/( Y[0,8] + X5 )**2
        Y[5,13]  = - D/( Y[0,8] + X5 )**2
        Y[6,13]  = - D*Y[6,8]/( Y[0,8] + X5 )**2  + 2.0*D*Y[1,8]**2/( Y[0,8] + X5 )**3
        Y[7,13]  = - D*Y[7,8]/( Y[0,8] + X5 )**2  + 2.0*D*Y[1,8]*Y[2,8]/( Y[0,8] + X5 )**3
        Y[8,13]  = - D*Y[8,8]/( Y[0,8] + X5 )**2  + 2.0*D*Y[2,8]**2/( Y[0,8] + X5 )**3
        Y[9,13]  = - D*Y[9,8]/( Y[0,8] + X5 )**2  + 2.0*D*Y[1,8]*Y[3,8]/( Y[0,8] + X5 )**3
        Y[10,13] = - D*Y[10,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[2,8]*Y[3,8]/( Y[0,8] + X5 )**3
        Y[11,13] = - D*Y[11,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[3,8]**2/( Y[0,8] + X5 )**3
        Y[12,13] = - D*Y[12,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[1,8]*Y[4,8]/( Y[0,8] + X5 )**3
        Y[13,13] = - D*Y[13,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[2,8]*Y[4,8]/( Y[0,8] + X5 )**3
        Y[14,13] = - D*Y[14,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[3,8]*Y[4,8]/( Y[0,8] + X5 )**3
        Y[15,13] = - D*Y[15,8]/( Y[0,8] + X5 )**2 + 2.0*D*Y[4,8]**2/( Y[0,8] + X5 )**3
        Y[16,13] = 2.0*D*Y[1,8]/( Y[0,8] + X5 )**3
        Y[17,13] = 2.0*D*Y[2,8]/( Y[0,8] + X5 )**3
        Y[18,13] = 2.0*D*Y[3,8]/( Y[0,8] + X5 )**3
        Y[19,13] = 2.0*D*Y[4,8]/( Y[0,8] + X5 )**3
        Y[20,13] = 2.0*D/( Y[0,8] + X5 )**3

        # Function C13.

        A        = 9.995e-1
        B        = 1.121e-1
        D        = 4.8e+1
        E        = 5.095e+3
        F        = 6.08e+1
        C[0,12]  = A*Y[0,9] - B*Y[0,13] + F*X2 + D*X4 - E
        C[1,12]  = A*Y[1,9] - B*Y[1,13]
        C[2,12]  = F + A*Y[2,9] - B*Y[2,13]
        C[3,12]  = - B*Y[3,13]
        C[4,12]  = D + A*Y[4,9] - B*Y[4,13]
        C[5,12]  = - B*Y[5,13]
        C[6,12]  = A*Y[6,9] - B*Y[6,13]
        C[7,12]  = A*Y[7,9] - B*Y[7,13]
        C[8,12]  = A*Y[8,9] - B*Y[8,13]
        C[9,12]  = - B*Y[9,13]
        C[10,12] = - B*Y[10,13]
        C[11,12] = - B*Y[11,13]
        C[12,12] = A*Y[12,9] - B*Y[12,13]
        C[13,12] = A*Y[13,9] - B*Y[13,13]
        C[14,12] = - B*Y[14,13]
        C[15,12] = A*Y[15,9] - B*Y[15,13]
        C[16,12] = - B*Y[16,13]
        C[17,12] = - B*Y[17,13]
        C[18,12] = - B*Y[18,13]
        C[19,12] = - B*Y[19,13]
        C[20,12] = - B*Y[20,13]

        # Function Y15.

        Y[0,14]  = Y[0,12]/C[0,12]
        Y[1,14]  = Y[1,12]/C[0,12] - Y[0,12]*C[1,12]/C[0,12]**2
        Y[2,14]  = Y[2,12]/C[0,12] - Y[0,12]*C[2,12]/C[0,12]**2
        Y[3,14]  = - Y[0,12]*C[3,12]/C[0,12]**2
        Y[4,14]  = Y[4,12]/C[0,12] - Y[0,12]*C[4,12]/C[0,12]**2
        Y[5,14]  = - Y[0,12]*C[5,12]/C[0,12]**2
        Y[6,14]  = ( Y[6,12]/C[0,12] - 2.0*Y[1,12]*C[1,12]/C[0,12]**2 
                   - Y[0,12]*C[6,12]/C[0,12]**2 + 2.0*Y[0,12]*C[1,12]**2/C[0,12]**3 )
        Y[7,14]  = ( Y[7,12]/C[0,12] - Y[1,12]*C[2,12]/C[0,12]**2 - Y[2,12]*C[1,12]/C[0,12]**2
                   - Y[0,12]*C[7,12]/C[0,12]**2 + 2.0*Y[0,12]*C[1,12]*C[2,12]/C[0,12]**3 )
        Y[8,14]  = ( Y[8,12]/C[0,12] - 2.0*Y[2,12]*C[2,12]/C[0,12]**2
                   - Y[0,12]*C[8,12]/C[0,12]**2 + 2.0*Y[0,12]*C[2,12]**2/C[0,12]**3 )
        Y[9,14]  = ( Y[9,12]/C[0,12] - Y[1,12]*C[3,12]/C[0,12]**2 - Y[3,12]*C[1,12]/C[0,12]**2
                   - Y[0,12]*C[9,12]/C[0,12]**2 + 2.0*Y[0,12]*C[1,12]*C[3,12]/C[0,12]**3 )
        Y[10,14] = ( Y[2,12]*C[3,12]/C[0,12]**2 - Y[0,12]*C[10,12]/C[0,12]**2 
                   + 2.0*Y[0,12]*C[2,12]*C[3,12]/C[0,12]**3 )
        Y[11,14] = - Y[0,12]*C[11,12]/C[0,12]**2 + 2.0*Y[0,12]*C[3,12]**2/C[0,12]**3
        Y[12,14] = ( Y[12,12]/C[0,12] - Y[1,12]*C[4,12]/C[0,12]**2 - Y[4,12]*C[1,12]/C[0,12]**2
                   - Y[0,12]*C[12,12]/C[0,12]**2 + 2.0*Y[0,12]*C[1,12]*C[4,12]/C[0,12]**3 )
        Y[13,14] = ( Y[13,12]/C[0,12] - Y[2,12]*C[4,12]/C[0,12]**2 - Y[4,12]*C[2,12]/C[0,12]**2
                   - Y[0,12]*C[13,12]/C[0,12]**2 + 2.0*Y[0,12]*C[2,12]*C[4,12]/C[0,12]**3 )
        Y[14,14] = ( - Y[4,12]*C[3,12]/C[0,12]**2 - Y[0,12]*C[14,12]/C[0,12]**2
                   + 2.0*Y[0,12]*C[3,12]*C[4,12]/C[0,12]**3 )
        Y[15,14] = ( Y[15,12]/C[0,12] - 2.0*Y[4,12]*C[4,12]/C[0,12]**2 
                   - Y[0,12]*C[15,12]/C[0,12]**2 + 2.0*Y[0,12]*C[4,12]**2/C[0,12]**3 )
        Y[16,14] = ( - Y[1,12]*C[5,12]/C[0,12]**2 -Y[0,12]*C[16,12]/C[0,12]**2
                   + 2.0*Y[0,12]*C[1,12]*C[5,12]/C[0,12]**3 )
        Y[17,14] = (- Y[2,12]*C[5,12]/C[0,12]**2 - Y[0,12]*C[17,12]/C[0,12]**2
                   + 2.0*Y[0,12]*C[2,12]*C[5,12]/C[0,12]**3 )
        Y[18,14] = - Y[0,12]*C[18,12]/C[0,12]**2 + 2.0*Y[0,12]*C[3,12]*C[5,12]/C[0,12]**3
        Y[19,14] = ( - Y[4,12]*C[5,12]/C[0,12]**2 - Y[0,12]*C[19,12]/C[0,12]**2
                   + 2.0*Y[0,12]*C[4,12]*C[5,12]/C[0,12]**3 )
        Y[20,14] = - Y[0,12]*C[20,12]/C[0,12]**2 + 2.0*Y[0,12]*C[5,12]**2/C[0,12]**3

        # Function Y16.

        A        = 1.48e+5
        B        = -3.31e+5
        D        = 4.0e+1
        E        = -6.1e+1
        Y[0,15]  = A + B*Y[0,14] + D*Y[0,12] + E*Y[0,14]*Y[0,12]
        Y[1,15]  = B*Y[1,14] + D*Y[1,12] + E*( Y[1,14]*Y[0,12] + Y[0,14]*Y[1,12])
        Y[2,15]  = B*Y[2,14] + D*Y[2,12] + E*( Y[2,14]*Y[0,12] + Y[0,14]*Y[2,12])
        Y[3,15]  = B*Y[3,14] + E*Y[3,14]*Y[0,12]
        Y[4,15]  = B*Y[4,14] + D*Y[4,12] + E*( Y[4,14]*Y[0,12] + Y[0,14]*Y[4,12])
        Y[5,15]  = B*Y[5,14] + E*Y[5,14]*Y[0,12]
        Y[6,15]  = ( B*Y[6,14] + D*Y[6,12]
                   + E*( Y[6,14]*Y[0,12] + Y[0,14]*Y[6,12]  + 2.0*Y[1,14]*Y[1,12] ) )
        Y[7,15]  = ( B*Y[7,14] + D*Y[7,12]
                   + E*( Y[7,14]*Y[0,12] + Y[0,14]*Y[7,12] + Y[1,14]*Y[2,12] + Y[2,14]*Y[1,12] ) )
        Y[8,15]  = ( B*Y[8,14] + D*Y[8,12] 
                   + E*( Y[8,14]*Y[0,12] + Y[0,14]*Y[8,12] + 2.0*Y[2,14]*Y[2,12] ) )
        Y[9,15]  = B*Y[9,14] + E*( Y[9,14]*Y[0,12] + Y[3,14]*Y[1,12])
        Y[10,15] = B*Y[10,14] + E*( Y[10,14]*Y[0,12] + Y[3,14]*Y[2,12] )
        Y[11,15] = B*Y[11,14] + E*Y[11,14]*Y[0,12]
        Y[12,15] = ( B*Y[12,14] + D*Y[12,12]
                   + E*( Y[12,14]*Y[0,12] + Y[0,14]*Y[12,12] + Y[1,14]*Y[4,12] + Y[4,14]*Y[1,12] ) )
        Y[13,15] = ( B*Y[13,14] + D*Y[13,12] 
                   + E*( Y[13,14]*Y[0,12] + Y[0,14]*Y[13,12] + Y[2,14]*Y[4,12] + Y[4,14]*Y[2,12] ) )
        Y[14,15] = B*Y[14,14] + E*( Y[14,14]*Y[0,12] + Y[3,14]*Y[4,12])
        Y[15,15] = ( B*Y[15,14] + D*Y[15,12]
                   + E*( Y[15,14]*Y[0,12] + Y[0,14]*Y[15,12] + 2.0*Y[4,14]*Y[4,12] ) )
        Y[16,15] = B*Y[16,14] + E*( Y[16,14]*Y[0,12] + Y[5,14]*Y[1,12])
        Y[17,15] = B*Y[17,14] + E*( Y[17,14]*Y[0,12] + Y[5,14]*Y[2,12])
        Y[18,15] = B*Y[18,14] + E*Y[18,14]*Y[0,12]
        Y[19,15] = B*Y[19,14] + E*( Y[19,14]*Y[0,12] + Y[5,14]*Y[4,12] )
        Y[20,15] = B*Y[20,14] + E*Y[20,14]*Y[0,12]

        # Function C14.

        A        = 2.324e+3
        B        = 2.874e+7
        C[0,13]  = A*Y[0,9] - B*Y[0,1]
        C[1,13]  = A*Y[1,9]
        C[2,13]  = A*Y[2,9]
        C[4,13]  = A*Y[4,9] - B*Y[4,1]
        C[6,13]  = A*Y[6,9]
        C[7,13]  = A*Y[7,9]
        C[8,13]  = A*Y[8,9]
        C[12,13] = A*Y[12,9]
        C[13,13] = A*Y[13,9]
        C[15,13] = A*Y[15,9] - B*Y[15,1]

        # Function Y17.

        A        = 1.413e+7
        B        = -1.328e+3
        D        = -5.31e+2
        Y[0,16]  = A + B*Y[0,9] + D*Y[0,10] + C[0,13]/C[0,11]
        Y[1,16]  = B*Y[1,9] + D*Y[1,10] + C[1,13]/C[0,11] - C[0,13]*C[1,11]/C[0,11]**2
        Y[2,16]  = B*Y[2,9] + C[2,13]/C[0,11] - C[0,13]*C[2,11]/C[0,11]**2
        Y[4,16]  = B*Y[4,9] + D*Y[4,10] + C[4,13]/C[0,11] - C[0,13]*C[4,11]/C[0,11]**2
        Y[6,16]  = ( B*Y[6,9] + D*Y[6,10] + C[6,13]/C[0,11] - C[0,13]*C[6,11]/C[0,11]**2
                   - 2.0*C[1,13]*C[1,11]/C[0,11]**2 + 2.0*C[0,13]*C[1,11]**2/C[0,11]**3 )
        Y[7,16]  = ( B*Y[7,9] + C[7,13]/C[0,11] - C[1,13]*C[2,11]/C[0,11]**2
                   - C[2,13]*C[1,11]/C[0,11]**2 - C[0,13]*C[7,11]/C[0,11]**2
                   + 2.0*C[0,13]*C[1,11]*C[2,11]/C[0,11]**3 )
        Y[8,16]  = ( B*Y[8,9] + C[8,13]/C[0,11] - C[0,13]*C[8,11]/C[0,11]**2
                   - 2.0*C[2,13]*C[2,11]/C[0,11]**2 + 2.0*C[0,13]*C[2,11]**2/C[0,11]**3 )
        Y[12,16] = ( B*Y[12,9] + D*Y[12,10] + C[12,13]/C[0,11] - C[1,13]*C[4,11]/C[0,11]**2
                   - C[4,13]*C[1,11]/C[0,11]**2 - C[0,13]*C[12,11]/C[0,11]**2
                   + 2.0*C[0,13]*C[1,11]*C[4,11]/C[0,11]**3 )
        Y[13,16] = ( B*Y[13,9] + C[13,13]/C[0,11] - C[2,13]*C[4,11]/C[0,11]**2
                   - C[4,13]*C[2,11]/C[0,11]**2
                   - C[0,13]*C[13,11]/C[0,11]**2 + 2.0*C[0,13]*C[2,11]*C[4,11]/C[0,11]**3 )
        Y[15,16] = ( B*Y[15,9] + D*Y[15,10] + C[15,13]/C[0,11] - C[0,13]*C[15,11]/C[0,11]**2
                   - 2.0*C[4,13]*C[4,11]/C[0,11]**2 + 2.0*C[0,13]*C[4,11]**2/C[0,11]**3 )

        # Function C15.

        A        = 5.2e-1
        C[0,14]  = Y[0,12]/Y[0,14] - Y[0,12]/A
        C[1,14]  = Y[1,12]/Y[0,14] - Y[0,12]*Y[1,14]/Y[0,14]**2 - Y[1,12]/A
        C[2,14]  = Y[2,12]/Y[0,14] - Y[0,12]*Y[2,14]/Y[0,14]**2 - Y[2,12]/A
        C[3,14]  = - Y[0,12]*Y[3,14]/Y[0,14]**2
        C[4,14]  = Y[4,12]/Y[0,14] - Y[0,12]*Y[4,14]/Y[0,14]**2 - Y[4,12]/A
        C[5,14]  = - Y[0,12]*Y[5,14]/Y[0,14]**2
        C[6,14]  = ( Y[6,12]/Y[0,14] - 2.0*Y[1,12]*Y[1,14]/Y[0,14]**2 - Y[0,12]*Y[6,14]/Y[0,14]**2
                   + 2.0*Y[0,12]*Y[1,14]**2/Y[0,14]**3 - Y[6,12]/A )
        C[7,14]  = ( Y[7,12]/Y[0,14] - Y[1,12]*Y[2,14]/Y[0,14]**2 - Y[2,12]*Y[1,14]/Y[0,14]**2
                   - Y[0,12]*Y[7,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[1,14]*Y[2,14]/Y[0,14]**3
                   - Y[7,12]/A )
        C[8,14]  = ( Y[8,12]/Y[0,14] - 2.0*Y[2,12]*Y[2,14]/Y[0,14]**2 - Y[0,12]*Y[8,14]/Y[0,14]**2
                   + 2.0*Y[0,12]*Y[2,14]**2/Y[0,14]**3 - Y[8,12]/A )
        C[9,14] = ( Y[9,12]/Y[0,14] - Y[1,12]*Y[3,14]/Y[0,14]**2 - Y[3,12]*Y[1,14]/Y[0,14]**2
                   - Y[0,12]*Y[9,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[1,14]*Y[3,14]/Y[0,14]**3 )
        C[10,14] = ( Y[2,12]*Y[3,14]/Y[0,14]**2 - Y[0,12]*Y[10,14]/Y[0,14]**2
                   + 2.0*Y[0,12]*Y[2,14]*Y[3,14]/Y[0,14]**3 )
        C[11,14] = - Y[0,12]*Y[11,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[3,14]**2/Y[0,14]**3
        C[12,14] = ( Y[12,12]/Y[0,14] - Y[1,12]*Y[4,14]/Y[0,14]**2 - Y[4,12]*Y[1,14]/Y[0,14]**2
                   - Y[0,12]*Y[12,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[1,14]*Y[4,14]/Y[0,14]**3
                   - Y[12,12]/A )
        C[13,14] = ( Y[13,12]/Y[0,14] - Y[2,12]*Y[4,14]/Y[0,14]**2 - Y[4,12]*Y[2,14]/Y[0,14]**2
                   - Y[0,12]*Y[13,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[2,14]*Y[4,14]/Y[0,14]**3
                   - Y[13,12]/A )
        C[14,14] = ( - Y[4,12]*Y[3,14]/Y[0,14]**2 - Y[0,12]*Y[14,14]/Y[0,14]**2 
                   + 2.0*Y[0,12]*Y[3,14]*Y[4,14]/Y[0,14]**3 )
        C[15,14] = ( Y[15,12]/Y[0,14] - 2.0*Y[4,12]*Y[4,14]/Y[0,14]**2
                   - Y[0,12]*Y[15,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[4,14]**2/Y[0,14]**3
                   - Y[15,12]/A )
        C[16,14] = ( - Y[1,12]*Y[5,14]/Y[0,14]**2 - Y[0,12]*Y[16,14]/Y[0,14]**2
                   + 2.0*Y[0,12]*Y[1,14]*Y[5,14]/Y[0,14]**3 )
        C[17,14] = ( - Y[2,12]*Y[5,14]/Y[0,14]**2 - Y[0,12]*Y[17,14]/Y[0,14]**2 
                   + 2.0*Y[0,12]*Y[2,14]*Y[5,14]/Y[0,14]**3 )
        C[18,14] = - Y[0,12]*Y[18,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[3,14]*Y[5,14]/Y[0,14]**3
        C[19,14] = ( - Y[4,12]*Y[5,14]/Y[0,14]**2 - Y[0,12]*Y[19,14]/Y[0,14]**2
                   + 2.0*Y[0,12]*Y[4,14]*Y[5,14]/Y[0,14]**3 )
        C[20,14] = - Y[0,12]*Y[20,14]/Y[0,14]**2 + 2.0*Y[0,12]*Y[5,14]**2/Y[0,14]**3

        # Function C16.

        A        = 1.104
        B        = 7.2e-1
        C[0,15]  = A - B*Y[0,14]
        C[1,15]  = - B*Y[1,14]
        C[2,15]  = - B*Y[2,14]
        C[3,15]  = - B*Y[3,14]
        C[4,15]  = - B*Y[4,14]
        C[5,15]  = - B*Y[5,14]
        C[6,15]  = - B*Y[6,14]
        C[7,15]  = - B*Y[7,14]
        C[8,15]  = - B*Y[8,14]
        C[9,15]  = - B*Y[9,14]
        C[10,15] = - B*Y[10,14]
        C[11,15] = - B*Y[11,14]
        C[12,15] = - B*Y[12,14]
        C[13,15] = - B*Y[13,14]
        C[14,15] = - B*Y[14,14]
        C[15,15] = - B*Y[15,14]
        C[16,15] = - B*Y[16,14]
        C[17,15] = - B*Y[17,14]
        C[18,15] = - B*Y[18,14]
        C[19,15] = - B*Y[19,14]
        C[20,15] = - B*Y[20,14]

        # Function C17.

        C[0,16]  = Y[0,8] + X5
        C[1,16]  = Y[1,8]
        C[2,16]  = Y[2,8]
        C[3,16]  = Y[3,8]
        C[4,16]  = Y[4,8]
        C[5,16]  = 1.0
        C[6,16]  = Y[6,8]
        C[7,16]  = Y[7,8]
        C[8,16]  = Y[8,8]
        C[9,16]  = Y[9,8]
        C[10,16] = Y[10,8]
        C[11,16] = Y[11,8]
        C[12,16] = Y[12,8]
        C[13,16] = Y[13,8]
        C[14,16] = Y[14,8]
        C[15,16] = Y[15,8]

        # Function C18.

        C[0,17]  = C[0,14]/C[0,15]
        C[1,17]  = C[1,14]/C[0,15] - C[0,14]*C[1,15]/C[0,15]**2
        C[2,17]  = C[2,14]/C[0,15] - C[0,14]*C[2,15]/C[0,15]**2
        C[3,17]  = C[3,14]/C[0,15] - C[0,14]*C[3,15]/C[0,15]**2
        C[4,17]  = C[4,14]/C[0,15] - C[0,14]*C[4,15]/C[0,15]**2
        C[5,17]  = C[5,14]/C[0,15] - C[0,14]*C[5,15]/C[0,15]**2
        C[6,17]  = ( C[6,14]/C[0,15] - 2.0*C[1,14]*C[1,15]/C[0,15]**2
                   - C[0,14]*C[6,15]/C[0,15]**2 + 2.0*C[0,14]*C[1,15]**2/C[0,15]**3 )
        C[7,17]  = ( C[7,14]/C[0,15] - C[1,14]*C[2,15]/C[0,15]**2 - C[2,14]*C[1,15]/C[0,15]**2
                   - C[0,14]*C[7,15]/C[0,15]**2 + 2.0*C[0,14]*C[1,15]*C[2,15]/C[0,15]**3 )
        C[8,17]  = ( C[8,14]/C[0,15] - 2.0*C[2,14]*C[2,15]/C[0,15]**2
                   - C[0,14]*C[8,15]/C[0,15]**2 + 2.0*C[0,14]*C[2,15]**2/C[0,15]**3 )
        C[9,17] = ( C[9,14]/C[0,15] - C[1,14]*C[3,15]/C[0,15]**2 - C[3,14]*C[1,15]/C[0,15]**2
                   - C[0,14]*C[9,15]/C[0,15]**2 +  2.0*C[0,14]*C[1,15]*C[3,15]/C[0,15]**3 )
        C[10,17] = ( C[10,14]/C[0,15] - C[2,14]*C[3,15]/C[0,15]**2 - C[3,14]*C[2,15]/C[0,15]**2
                   - C[0,14]*C[10,15]/C[0,15]**2 + 2.0*C[0,14]*C[2,15]*C[3,15]/C[0,15]**3 )
        C[11,17] = ( C[11,14]/C[0,15] - 2.0*C[3,14]*C[3,15]/C[0,15]**2
                   - C[0,14]*C[11,15]/C[0,15]**2 + 2.0*C[0,14]*C[3,15]**2/C[0,15]**3 )
        C[12,17] = ( C[12,14]/C[0,15] - C[1,14]*C[4,15]/C[0,15]**2 - C[4,14]*C[1,15]/C[0,15]**2
                   - C[0,14]*C[12,15]/C[0,15]**2 + 2.0*C[0,14]*C[1,15]*C[4,15]/C[0,15]**3 )
        C[13,17] = ( C[13,14]/C[0,15] - C[2,14]*C[4,15]/C[0,15]**2 - C[4,14]*C[2,15]/C[0,15]**2
                   - C[0,14]*C[13,15]/C[0,15]**2 + 2.0*C[0,14]*C[2,15]*C[4,15]/C[0,15]**3 )
        C[14,17] = ( C[14,14]/C[0,15] - C[3,14]*C[4,15]/C[0,15]**2 - C[4,14]*C[3,15]/C[0,15]**2
                   - C[0,14]*C[14,15]/C[0,15]**2 + 2.0*C[0,14]*C[3,15]*C[4,15]/C[0,15]**3 )
        C[15,17] = ( C[15,14]/C[0,15] - 2.0*C[4,14]*C[4,15]/C[0,15]**2
                   - C[0,14]*C[15,15]/C[0,15]**2 + 2.0*C[0,14]*C[4,15]**2/C[0,15]**3 )
        C[16,17] = ( C[16,14]/C[0,15] - C[1,14]*C[5,15]/C[0,15]**2 - C[5,14]*C[1,15]/C[0,15]**2
                   - C[0,14]*C[16,15]/C[0,15]**2 + 2.0*C[0,14]*C[1,15]*C[5,15]/C[0,15]**3 )
        C[17,17] = ( C[17,14]/C[0,15] - C[2,14]*C[5,15]/C[0,15]**2 - C[5,14]*C[2,15]/C[0,15]**2
                   - C[0,14]*C[17,15]/C[0,15]**2 + 2.0*C[0,14]*C[2,15]*C[5,15]/C[0,15]**3 )
        C[18,17] = ( C[18,14]/C[0,15] - C[3,14]*C[5,15]/C[0,15]**2 - C[5,14]*C[3,15]/C[0,15]**2
                   - C[0,14]*C[18,15]/C[0,15]**2 + 2.0*C[0,14]*C[3,15]*C[5,15]/C[0,15]**3 )
        C[19,17] = ( C[19,14]/C[0,15] - C[4,14]*C[5,15]/C[0,15]**2 - C[5,14]*C[4,15]/C[0,15]**2
                   - C[0,14]*C[19,15]/C[0,15]**2 + 2.0*C[0,14]*C[4,15]*C[5,15]/C[0,15]**3 )
        C[20,17] = ( C[20,14]/C[0,15] - 2.0*C[5,14]*C[5,15]/C[0,15]**2
                   - C[0,14]*C[20,15]/C[0,15]**2 + 2.0*C[0,14]*C[5,15]**2/C[0,15]**3 )

        # Function C19.

        C[0,18]  = Y[0,1]/C[0,11]
        C[1,18]  = - Y[0,1]*C[1,11]/ C[0,11]**2
        C[2,18]  = - Y[0,1]*C[2,11]/ C[0,11]**2
        C[4,18]  = Y[4,1]/C[0,11] - Y[0,1]*C[4,11]/ C[0,11]**2
        C[6,18]  = Y[0,1]*( 2.0*C[1,11]**2/C[0,11]**3 - C[6,11]/ C[0,11]**2)
        C[7,18]  = Y[0,1]*( 2.0*C[1,11]*C[2,11]/C[0,11]**3 - C[7,11]/ C[0,11]**2 )
        C[8,18]  = Y[0,1]*( 2.0*C[2,11]**2/C[0,11]**3 - C[8,11]/ C[0,11]**2)
        C[12,18] = ( Y[0,1]*( 2.0*C[1,11]*C[4,11]/C[0,11]**3 - C[12,11]/C[0,11]**2 ) 
                   - Y[4,1]*C[1,11]/C[0,11]**2 )
        C[13,18] = ( Y[0,1]*( 2.0*C[2,11]*C[4,11]/C[0,11]**3 - C[13,11]/C[0,11]**2 ) 
                   - Y[4,1]*C[2,11]/C[0,11]**2 )
        C[15,18] = ( Y[15,1]/C[0,11] - 2.0*Y[4,1]*C[4,11]/C[0,11]**2
                   - Y[0,1]*C[15,11]/ C[0,11]**2 + 2.0*Y[0,1]*C[4,11]**2/ C[0,11]**3 )

        # Function Y18.

        A        = - 2.8e-1/7.2e-1
        Y[0,17]  = Y[0,3] + A/Y[0,4]
        Y[1,17]  = Y[1,3] - A*Y[1,4]/Y[0,4]**2
        Y[2,17]  = - A*Y[2,4]/Y[0,4]**2
        Y[4,17]  = Y[4,3] - A*Y[4,4]/Y[0,4]**2
        Y[6,17]  = Y[6,3] - A*Y[6,4]/Y[0,4]**2 + 2.0*A*Y[1,4]**2/Y[0,4]**3
        Y[7,17]  = - A*Y[7,4]/Y[0,4]**2 + 2.0*A*Y[1,4]*Y[2,4]/Y[0,4]**3
        Y[8,17]  = - A*Y[8,4]/Y[0,4]**2 + 2.0*A*Y[2,4]**2/Y[0,4]**3
        Y[12,17] = Y[12,3] - A*Y[12,4]/Y[0,4]**2 + 2.0*A*Y[1,4]*Y[4,4]/Y[0,4]**3
        Y[13,17] = - A*Y[13,4]/Y[0,4]**2 + 2.0*A*Y[2,4]*Y[4,4]/Y[0,4]**3
        Y[15,17] = Y[15,3] - A*Y[15,4]/Y[0,4]**2 + 2.0*A*Y[4,4]**2/Y[0,4]**3

        # Function Y19.

        A        = -3.496e+3
        Y[0,18]  = A*C[0,18]
        Y[1,18]  = A*C[1,18]
        Y[2,18]  = A*C[2,18]
        Y[4,18]  = A*C[4,18]
        Y[6,18]  = A*C[6,18]
        Y[7,18]  = A*C[7,18]
        Y[8,18]  = A*C[8,18]
        Y[12,18] = A*C[12,18]
        Y[13,18] = A*C[13,18]
        Y[15,18] = A*C[15,18]

        # Function Y20

        A        = 6.2212e+4
        Y[0,19]  =   A/C[0,16] - Y[0,0]
        Y[1,19]  = - A*C[1,16] /C[0,16]**2
        Y[2,19]  = - A*C[2,16] /C[0,16]**2 - Y[2,0]
        Y[3,19]  = - A*C[3,16] /C[0,16]**2 - Y[3,0]
        Y[4,19]  = - A*C[4,16] /C[0,16]**2
        Y[5,19]  = - A*C[5,16] /C[0,16]**2
        Y[6,19]  = - A*C[6,16] /C[0,16]**2 + 2.0*A*C[1,16]**2/C[0,16]**3
        Y[7,19]  = - A*C[7,16] /C[0,16]**2 + 2.0*A*C[1,16]*C[2,16]/C[0,16]**3
        Y[8,19]  = - A*C[8,16] /C[0,16]**2 + 2.0*A*C[2,16]**2/C[0,16]**3
        Y[9,19]  = - A*C[9,16] /C[0,16]**2 + 2.0*A*C[1,16]*C[3,16]/C[0,16]**3
        Y[10,19] = - A*C[10,16]/C[0,16]**2 + 2.0*A*C[2,16]*C[3,16]/C[0,16]**3
        Y[11,19] = - A*C[11,16]/C[0,16]**2 + 2.0*A*C[3,16]**2/C[0,16]**3
        Y[12,19] = - A*C[12,16]/C[0,16]**2 + 2.0*A*C[1,16]*C[4,16]/C[0,16]**3
        Y[13,19] = - A*C[13,16]/C[0,16]**2 + 2.0*A*C[2,16]*C[4,16]/C[0,16]**3
        Y[14,19] = - A*C[14,16]/C[0,16]**2 + 2.0*A*C[3,16]*C[4,16]/C[0,16]**3
        Y[15,19] = - A*C[15,16]/C[0,16]**2 + 2.0*A*C[4,16]**2/C[0,16]**3
        Y[16,19] = - A*C[16,16]/C[0,16]**2 + 2.0*A*C[1,16]*C[5,16]/C[0,16]**3
        Y[17,19] = - A*C[17,16]/C[0,16]**2 + 2.0*A*C[2,16]*C[5,16]/C[0,16]**3
        Y[18,19] = - A*C[18,16]/C[0,16]**2 + 2.0*A*C[3,16]*C[5,16]/C[0,16]**3
        Y[19,19] = - A*C[19,16]/C[0,16]**2 + 2.0*A*C[4,16]*C[5,16]/C[0,16]**3
        Y[20,19] = - A*C[20,16]/C[0,16]**2 + 2.0*A*C[5,16]**2/C[0,16]**3
      
        return C, Y
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

