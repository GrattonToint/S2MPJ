from s2mpjlib import *
class  CSFI2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CSFI2
#    *********
# 
#    Source: problem MINLEN in
#    Vasko and Stott
#    "Optimizing continuous caster product dimensions:
#     an example of a nonlinear design problem in the steel industry"
#    SIAM Review, Vol 37 No, 1 pp.82-84, 1995
# 
#    SIF input: A.R. Conn April 1995
# 
#    classification = "C-CLOR2-RN-5-4"
# 
#    input parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CSFI2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['MINTPH'] = 45.0
        v_['MINTHICK'] = 7.0
        v_['MINAREA'] = 200.0
        v_['MAXAREA'] = 250.0
        v_['MAXASPR'] = 2.0
        v_['K'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('THICK',ix_)
        self.xnames=arrset(self.xnames,iv,'THICK')
        [iv,ix_,_] = s2mpj_ii('WID',ix_)
        self.xnames=arrset(self.xnames,iv,'WID')
        [iv,ix_,_] = s2mpj_ii('LEN',ix_)
        self.xnames=arrset(self.xnames,iv,'LEN')
        [iv,ix_,_] = s2mpj_ii('TPH',ix_)
        self.xnames=arrset(self.xnames,iv,'TPH')
        [iv,ix_,_] = s2mpj_ii('IPM',ix_)
        self.xnames=arrset(self.xnames,iv,'IPM')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['LEN']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CIPM',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CIPM')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['IPM']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('CLEN',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CLEN')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['LEN']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('WOT',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'WOT')
        [ig,ig_,_] = s2mpj_ii('TTW',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'TTW')
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
        self.gconst = arrset(self.gconst,ig_['WOT'],float(v_['MAXASPR']))
        self.gconst = arrset(self.gconst,ig_['TTW'],float(v_['MINAREA']))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = -np.full((self.nle,1),float('inf'))
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        v_['RHS'] = v_['MAXAREA']-v_['MINAREA']
        grange = arrset(grange,ig_['TTW'],float(v_['RHS']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['THICK']] = v_['MINTHICK']
        self.xlower[ix_['TPH']] = v_['MINTPH']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eCMPLQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        [it,iet_,_] = s2mpj_ii( 'eSQQUT', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'eQUOTN', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCMPLQ')
        ielftype = arrset(ielftype,ie,iet_["eCMPLQ"])
        vname = 'TPH'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'WID'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'THICK'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQUT')
        ielftype = arrset(ielftype,ie,iet_["eSQQUT"])
        vname = 'THICK'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'IPM'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQUOTN')
        ielftype = arrset(ielftype,ie,iet_["eQUOTN"])
        vname = 'WID'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'THICK'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'THICK'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'WID'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['CIPM']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CLEN']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['WOT']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['TTW']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               55.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle)] = grange[legrps]
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nle+self.neq,self.m)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-RN-5-4"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCMPLQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP0 = EV_[1,0]*EV_[2,0]
        TMP1 = 117.3708920187793427e0*EV_[0,0]/TMP0
        TMP2 = 117.3708920187793427e0/TMP0
        f_   = TMP1
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TMP2
            g_[1] = -TMP1/EV_[1,0]
            g_[2] = -TMP1/EV_[2,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = -TMP2/EV_[1,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = -TMP2/EV_[2,0]
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0e0*TMP1/(EV_[1,0]*EV_[1,0])
                H_[1,2] = TMP1/TMP0
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0e0*TMP1/(EV_[2,0]*EV_[2,0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQQUT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = EV_[0,0]*EV_[1,0]/48.0e0
        f_   = EV_[0,0]*TMP
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*TMP
            g_[1] = EV_[0,0]*EV_[0,0]/48.0e0
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EV_[1,0]/24.0e0
                H_[0,1] = EV_[0,0]/24.0e0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQUOTN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = EV_[0,0]/EV_[1,0]
        f_   = TMP
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0e0/EV_[1,0]
            g_[1] = -TMP/EV_[1,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0e0/(EV_[1,0]*EV_[1,0])
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0e0*TMP/(EV_[1,0]*EV_[1,0])
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
                H_[0,1] = 1.0e0
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

