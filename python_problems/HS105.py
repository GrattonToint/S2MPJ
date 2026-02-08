from s2mpjlib import *
class  HS105(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS105
#    *********
# 
#    Source: problem 105 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
#    bug correction (line 351) Ph. Toint, May 2024
# 
#    classification = "C-COLR2-AY-8-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS105'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 8
        v_['1'] = 1
        v_['235'] = 235
        v_['Y1'] = 95.0
        v_['Y2'] = 105.0
        v_['LOW'] = 3
        v_['UP'] = 6
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 110.0
        v_['LOW'] = 7
        v_['UP'] = 10
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 115.0
        v_['LOW'] = 11
        v_['UP'] = 25
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 120.0
        v_['LOW'] = 26
        v_['UP'] = 40
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 125.0
        v_['LOW'] = 41
        v_['UP'] = 55
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 130.0
        v_['LOW'] = 56
        v_['UP'] = 68
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 135.0
        v_['LOW'] = 69
        v_['UP'] = 89
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 140.0
        v_['LOW'] = 90
        v_['UP'] = 101
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 145.0
        v_['LOW'] = 102
        v_['UP'] = 118
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 150.0
        v_['LOW'] = 119
        v_['UP'] = 122
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 155.0
        v_['LOW'] = 123
        v_['UP'] = 142
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 160.0
        v_['LOW'] = 143
        v_['UP'] = 150
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 165.0
        v_['LOW'] = 151
        v_['UP'] = 167
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 170.0
        v_['LOW'] = 168
        v_['UP'] = 175
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 175.0
        v_['LOW'] = 176
        v_['UP'] = 181
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 180.0
        v_['LOW'] = 182
        v_['UP'] = 187
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 185.0
        v_['LOW'] = 188
        v_['UP'] = 194
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 190.0
        v_['LOW'] = 195
        v_['UP'] = 198
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 195.0
        v_['LOW'] = 199
        v_['UP'] = 201
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 200.0
        v_['LOW'] = 202
        v_['UP'] = 204
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 205.0
        v_['LOW'] = 205
        v_['UP'] = 212
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 210.0
        v_['Y213'] = 215.0
        v_['LOW'] = 214
        v_['UP'] = 219
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 220.0
        v_['LOW'] = 220
        v_['UP'] = 224
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 230.0
        v_['Y225'] = 235.0
        v_['LOW'] = 226
        v_['UP'] = 232
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 240.0
        v_['Y233'] = 245.0
        v_['LOW'] = 234
        v_['UP'] = 235
        for I in range(int(v_['LOW']),int(v_['UP'])+1):
            v_['Y'+str(I)] = 250.0
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
        for I in range(int(v_['1']),int(v_['235'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0e+0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0e+0))
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
        self.gconst = arrset(self.gconst,ig_['C1'],float(-1.0e+0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X1']] = 0.001
        self.xupper[ix_['X1']] = 0.499
        self.xlower[ix_['X2']] = 0.001
        self.xupper[ix_['X2']] = 0.499
        self.xlower[ix_['X3']] = 100.0
        self.xupper[ix_['X3']] = 180.0
        self.xlower[ix_['X4']] = 130.0
        self.xupper[ix_['X4']] = 210.0
        self.xlower[ix_['X5']] = 170.0
        self.xupper[ix_['X5']] = 240.0
        self.xlower[ix_['X6']] = 5.0
        self.xupper[ix_['X6']] = 25.0
        self.xlower[ix_['X7']] = 5.0
        self.xupper[ix_['X7']] = 25.0
        self.xlower[ix_['X8']] = 5.0
        self.xupper[ix_['X8']] = 25.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(0.1)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(0.1)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(0.2)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(0.2))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(100.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(100.0)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(125.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(125.0))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(175.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(175.0)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(11.2)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X6'])[0],float(11.2))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(13.2)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X7']),float(13.2)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(15.8)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X8'])[0],float(15.8))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eABI', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'YI')
        [it,iet_,_] = s2mpj_ii( 'eCI', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X5')
        elftv = loaset(elftv,it,3,'X8')
        elftp = loaset(elftp,it,0,'YI')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['235'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eABI')
            ielftype = arrset(ielftype,ie,iet_["eABI"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='YI')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eABI')
            ielftype = arrset(ielftype,ie,iet_["eABI"])
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='YI')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eCI')
            ielftype = arrset(ielftype,ie,iet_["eCI"])
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X8'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='YI')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['235'])+1):
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gLOG')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
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
        self.pbclass   = "C-COLR2-AY-8-1"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eABI(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R = EV_[0,0]/EV_[1,0]
        D = (self.elpar[iel_][0]-EV_[2,0])/EV_[1,0]
        E = np.exp(-5.0e-1*D*D)
        DDV2 = -D/EV_[1,0]
        DEV2 = E*(-D)*DDV2
        DDV3 = -1.0e+0/EV_[1,0]
        DEV3 = E*(-D)*DDV3
        f_   = R*E
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/EV_[1,0]
            g_[1] = (D*D-1.0e+0)*R*E/EV_[1,0]
            g_[2] = D*R*E/EV_[1,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEV2-E/EV_[1,0])/EV_[1,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = DEV3/EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/EV_[1,0]))*R/EV_[1,0]
                H_[1,2] = (DDV2*E+D*DEV2-2.0e+0*D*E/EV_[1,0])*R/EV_[1,0]
                H_[2,1] = H_[1,2]
                H_[2,2] = (DDV3*E+D*DEV3)*R/EV_[1,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCI(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]-1
        U_[0,1] = U_[0,1]-1
        U_[1,3] = U_[1,3]+1
        U_[2,2] = U_[2,2]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        R = (1.0e+0+IV_[0])/IV_[1]
        D = (self.elpar[iel_][0]-IV_[2])/IV_[1]
        E = np.exp(-5.0e-1*D*D)
        DDV2 = -D/IV_[1]
        DEV2 = E*(-D)*DDV2
        DDV3 = -1.0e+0/IV_[1]
        DEV3 = E*(-D)*DDV3
        f_   = R*E
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E/IV_[1]
            g_[1] = (D*D-1.0e+0)*R*E/IV_[1]
            g_[2] = D*R*E/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (DEV2-E/IV_[1])/IV_[1]
                H_[1,0] = H_[0,1]
                H_[0,2] = DEV3/IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0e+0*D*DDV2*E+(D*D-1.0e+0)*(DEV2-2.0e+0*E/IV_[1]))*R/IV_[1]
                H_[1,2] = (DDV2*E+D*DEV2-2.0e+0*D*E/IV_[1])*R/IV_[1]
                H_[2,1] = H_[1,2]
                H_[2,2] = (DDV3*E+D*DEV3)*R/IV_[1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def g_globs(self):

        self.gfpar = np.array([]);
        self.gfpar = arrset(self.gfpar,0,3.9894228040143270e-01)     # this is P
        return pbm

    @staticmethod
    def gLOG(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= -np.log(self.gfpar[0]*GVAR_)
        if nargout>1:
            g_ = -1/GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 1/GVAR_**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

