from s2mpjlib import *
class  HS99(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS99
#    *********
# 
#    Source: problem 99 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-COOR2-AN-7-2"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS99'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['2'] = 2
        v_['7'] = 7
        v_['8'] = 8
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['7'])+1):
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
        self.gscale = arrset(self.gscale,ig,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('Q8E',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'Q8E')
        [ig,ig_,_] = s2mpj_ii('S8E',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'S8E')
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
        self.gconst = arrset(self.gconst,ig_['Q8E'],float(100000.0))
        self.gconst = arrset(self.gconst,ig_['S8E'],float(1000.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),1.58)
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eR8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        [it,iet_,_] = s2mpj_ii( 'eQ8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        [it,iet_,_] = s2mpj_ii( 'eS8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'R8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eR8T')
        ielftype = arrset(ielftype,ie,iet_["eR8T"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'Q8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQ8T')
        ielftype = arrset(ielftype,ie,iet_["eQ8T"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eS8T')
        ielftype = arrset(ielftype,ie,iet_["eS8T"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.58),float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='X7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        self.grftype = arrset(self.grftype,ig,'gL2')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['R8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['Q8E']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['Q8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['S8E']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -831079892.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AN-7-2"
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
        self.efpar = arrset( self.efpar,0,50.0)
        self.efpar = arrset( self.efpar,1,50.0)
        self.efpar = arrset( self.efpar,2,75.0)
        self.efpar = arrset( self.efpar,3,75.0)
        self.efpar = arrset( self.efpar,4,75.0)
        self.efpar = arrset( self.efpar,5,100.0)
        self.efpar = arrset( self.efpar,6,100.0)
        self.efpar = arrset( self.efpar,7,25.0)
        self.efpar = arrset( self.efpar,8,25.0)
        self.efpar = arrset( self.efpar,9,50.0)
        self.efpar = arrset( self.efpar,10,50.0)
        self.efpar = arrset( self.efpar,11,50.0)
        self.efpar = arrset( self.efpar,12,90.0)
        self.efpar = arrset( self.efpar,13,90.0)
        self.efpar = arrset( self.efpar,14,32.0)
        return pbm

    @staticmethod
    def eR8T(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R2 = self.efpar[0]*self.efpar[7]*np.cos(EV_[0,0])
        R3 = self.efpar[1]*self.efpar[8]*np.cos(EV_[1,0])+R2
        R4 = self.efpar[2]*self.efpar[9]*np.cos(EV_[2,0])+R3
        R5 = self.efpar[3]*self.efpar[10]*np.cos(EV_[3,0])+R4
        R6 = self.efpar[4]*self.efpar[11]*np.cos(EV_[4,0])+R5
        R7 = self.efpar[5]*self.efpar[12]*np.cos(EV_[5,0])+R6
        f_   = self.efpar[6]*self.efpar[13]*np.cos(EV_[6,0])+R7
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -self.efpar[0]*self.efpar[7]*np.sin(EV_[0,0])
            g_[1] = -self.efpar[1]*self.efpar[8]*np.sin(EV_[1,0])
            g_[2] = -self.efpar[2]*self.efpar[9]*np.sin(EV_[2,0])
            g_[3] = -self.efpar[3]*self.efpar[10]*np.sin(EV_[3,0])
            g_[4] = -self.efpar[4]*self.efpar[11]*np.sin(EV_[4,0])
            g_[5] = -self.efpar[5]*self.efpar[12]*np.sin(EV_[5,0])
            g_[6] = -self.efpar[6]*self.efpar[13]*np.sin(EV_[6,0])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = -self.efpar[0]*self.efpar[7]*np.cos(EV_[0,0])
                H_[1,1] = -self.efpar[1]*self.efpar[8]*np.cos(EV_[1,0])
                H_[2,2] = -self.efpar[2]*self.efpar[9]*np.cos(EV_[2,0])
                H_[3,3] = -self.efpar[3]*self.efpar[10]*np.cos(EV_[3,0])
                H_[4,4] = -self.efpar[4]*self.efpar[11]*np.cos(EV_[4,0])
                H_[5,5] = -self.efpar[5]*self.efpar[12]*np.cos(EV_[5,0])
                H_[6,6] = -self.efpar[6]*self.efpar[13]*np.cos(EV_[6,0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eS8T(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S2 = self.efpar[7]*(self.efpar[0]*np.sin(EV_[0,0])-self.efpar[14])
        S3 = self.efpar[8]*(self.efpar[1]*np.sin(EV_[1,0])-self.efpar[14])+S2
        S4 = self.efpar[9]*(self.efpar[2]*np.sin(EV_[2,0])-self.efpar[14])+S3
        S5 = self.efpar[10]*(self.efpar[3]*np.sin(EV_[3,0])-self.efpar[14])+S4
        S6 = self.efpar[11]*(self.efpar[4]*np.sin(EV_[4,0])-self.efpar[14])+S5
        S7 = self.efpar[12]*(self.efpar[5]*np.sin(EV_[5,0])-self.efpar[14])+S6
        f_   = self.efpar[13]*(self.efpar[6]*np.sin(EV_[6,0])-self.efpar[14])+S7
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.efpar[0]*self.efpar[7]*np.cos(EV_[0,0])
            g_[1] = self.efpar[1]*self.efpar[8]*np.cos(EV_[1,0])
            g_[2] = self.efpar[2]*self.efpar[9]*np.cos(EV_[2,0])
            g_[3] = self.efpar[3]*self.efpar[10]*np.cos(EV_[3,0])
            g_[4] = self.efpar[4]*self.efpar[11]*np.cos(EV_[4,0])
            g_[5] = self.efpar[5]*self.efpar[12]*np.cos(EV_[5,0])
            g_[6] = self.efpar[6]*self.efpar[13]*np.cos(EV_[6,0])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = -self.efpar[0]*self.efpar[7]*np.sin(EV_[0,0])
                H_[1,1] = -self.efpar[1]*self.efpar[8]*np.sin(EV_[1,0])
                H_[2,2] = -self.efpar[2]*self.efpar[9]*np.sin(EV_[2,0])
                H_[3,3] = -self.efpar[3]*self.efpar[10]*np.sin(EV_[3,0])
                H_[4,4] = -self.efpar[4]*self.efpar[11]*np.sin(EV_[4,0])
                H_[5,5] = -self.efpar[5]*self.efpar[12]*np.sin(EV_[5,0])
                H_[6,6] = -self.efpar[6]*self.efpar[13]*np.sin(EV_[6,0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQ8T(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S2 = self.efpar[7]*(self.efpar[0]*np.sin(EV_[0,0])-self.efpar[14])
        S3 = self.efpar[8]*(self.efpar[1]*np.sin(EV_[1,0])-self.efpar[14])+S2
        S4 = self.efpar[9]*(self.efpar[2]*np.sin(EV_[2,0])-self.efpar[14])+S3
        S5 = self.efpar[10]*(self.efpar[3]*np.sin(EV_[3,0])-self.efpar[14])+S4
        S6 = self.efpar[11]*(self.efpar[4]*np.sin(EV_[4,0])-self.efpar[14])+S5
        S7 = self.efpar[12]*(self.efpar[5]*np.sin(EV_[5,0])-self.efpar[14])+S6
        DSD1 = self.efpar[0]*self.efpar[7]*np.cos(EV_[0,0])
        DSD2 = self.efpar[1]*self.efpar[8]*np.cos(EV_[1,0])
        DSD3 = self.efpar[2]*self.efpar[9]*np.cos(EV_[2,0])
        DSD4 = self.efpar[3]*self.efpar[10]*np.cos(EV_[3,0])
        DSD5 = self.efpar[4]*self.efpar[11]*np.cos(EV_[4,0])
        DSD6 = self.efpar[5]*self.efpar[12]*np.cos(EV_[5,0])
        DSD7 = self.efpar[6]*self.efpar[13]*np.cos(EV_[6,0])
        D2SD1 = -self.efpar[0]*self.efpar[7]*np.sin(EV_[0,0])
        D2SD2 = -self.efpar[1]*self.efpar[8]*np.sin(EV_[1,0])
        D2SD3 = -self.efpar[2]*self.efpar[9]*np.sin(EV_[2,0])
        D2SD4 = -self.efpar[3]*self.efpar[10]*np.sin(EV_[3,0])
        D2SD5 = -self.efpar[4]*self.efpar[11]*np.sin(EV_[4,0])
        D2SD6 = -self.efpar[5]*self.efpar[12]*np.sin(EV_[5,0])
        D2SD7 = -self.efpar[6]*self.efpar[13]*np.sin(EV_[6,0])
        Q2  = (
              0.5*self.efpar[7]*self.efpar[7]*(self.efpar[0]*np.sin(EV_[0,0])-self.efpar[14]))
        Q3  = (
              0.5*self.efpar[8]*self.efpar[8]*(self.efpar[1]*np.sin(EV_[1,0])-self.efpar[14])+self.efpar[8]*S2+Q2)
        Q4  = (
              0.5*self.efpar[9]*self.efpar[9]*(self.efpar[2]*np.sin(EV_[2,0])-self.efpar[14])+self.efpar[9]*S3+Q3)
        Q5  = (
              0.5*self.efpar[10]*self.efpar[10]*(self.efpar[3]*np.sin(EV_[3,0])-self.efpar[14])+self.efpar[10]*S4+Q4)
        Q6  = (
              0.5*self.efpar[11]*self.efpar[11]*(self.efpar[4]*np.sin(EV_[4,0])-self.efpar[14])+self.efpar[11]*S5+Q5)
        Q7  = (
              0.5*self.efpar[12]*self.efpar[12]*(self.efpar[5]*np.sin(EV_[5,0])-self.efpar[14])+self.efpar[12]*S6+Q6)
        f_    = (
              0.5*self.efpar[13]*self.efpar[13]*(self.efpar[6]*np.sin(EV_[6,0])-self.efpar[14])+self.efpar[13]*S7+Q7)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (0.5*self.efpar[7]*self.efpar[7]*self.efpar[0]*np.cos(EV_[0,0])+
                 (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10]+self.efpar[9]+self.efpar[8])*DSD1)
            g_[1] = (0.5*self.efpar[8]*self.efpar[8]*self.efpar[1]*np.cos(EV_[1,0])+
                 (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10]+self.efpar[9])*DSD2)
            g_[2] = (0.5*self.efpar[9]*self.efpar[9]*self.efpar[2]*np.cos(EV_[2,0])+
                 (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10])*DSD3)
            g_[3] = (0.5*self.efpar[10]*self.efpar[10]*self.efpar[3]*np.cos(EV_[3,0])+
                 (self.efpar[13]+self.efpar[12]+self.efpar[11])*DSD4)
            g_[4] = (0.5*self.efpar[11]*self.efpar[11]*self.efpar[4]*np.cos(EV_[4,0])+
                 (self.efpar[13]+self.efpar[12])*DSD5)
            g_[5] = (0.5*self.efpar[12]*self.efpar[12]*self.efpar[5]*np.cos(EV_[5,0])+
                 self.efpar[13]*DSD6)
            g_[6] = 0.5*self.efpar[13]*self.efpar[13]*self.efpar[6]*np.cos(EV_[6,0])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = (-0.5*self.efpar[7]*self.efpar[7]*self.efpar[0]*np.sin(EV_[0,0])+
                     (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10]+self.efpar[9]+self.efpar[8])*D2SD1)
                H_[1,1] = (-0.5*self.efpar[8]*self.efpar[8]*self.efpar[1]*np.sin(EV_[1,0])+
                     (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10]+self.efpar[9])*D2SD2)
                H_[2,2] = (-0.5*self.efpar[9]*self.efpar[9]*self.efpar[2]*np.sin(EV_[2,0])+
                     (self.efpar[13]+self.efpar[12]+self.efpar[11]+self.efpar[10])*D2SD3)
                H_[3,3] = (-0.5*self.efpar[10]*self.efpar[10]*self.efpar[3]*np.sin(EV_[3,0])+
                     (self.efpar[13]+self.efpar[12]+self.efpar[11])*D2SD4)
                H_[4,4] = (-0.5*self.efpar[11]*self.efpar[11]*self.efpar[4]*np.sin(EV_[4,0])+
                     (self.efpar[13]+self.efpar[12])*D2SD5)
                H_[5,5] = (-0.5*self.efpar[12]*self.efpar[12]*self.efpar[5]*np.sin(EV_[5,0])+
                     self.efpar[13]*D2SD6)
                H_[6,6] = -0.5*self.efpar[13]*self.efpar[13]*self.efpar[6]*np.sin(EV_[6,0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

