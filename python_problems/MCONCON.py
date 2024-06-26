from s2mpjlib import *
class  MCONCON(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Another small gas network problem.
# 
#    SIF input: Sybille Schachler, Oxford, August 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "LOI2-MN-15-11"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MCONCON'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 7
        v_['M'] = 4
        v_['M+1'] = 1+v_['M']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'P'+str(I))
            [iv,ix_,_] = s2mpj_ii('Q'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Q'+str(I))
            [iv,ix_,_] = s2mpj_ii('F'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'F'+str(I))
        for I in range(int(v_['M+1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'P'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJECT',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['P'+str(I)]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('PAN'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PAN'+str(I))
        [ig,ig_,_] = s2mpj_ii('MBAL1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL1')
        iv = ix_['Q1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['F3']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL2')
        iv = ix_['Q1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['F1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL3')
        iv = ix_['Q2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['F1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL4')
        iv = ix_['Q2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['Q3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL5')
        iv = ix_['Q3']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['F2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL6')
        iv = ix_['Q4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['F2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('MBAL7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'MBAL7')
        iv = ix_['Q4']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['F4']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['DEMAND'] = -1000.0
        self.gconst = arrset(self.gconst,ig_['MBAL4'],float(v_['DEMAND']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        v_['PMAX1'] = 914.73
        v_['PMAX2'] = 904.73
        self.xupper[ix_['P3']] = v_['PMAX2']
        self.xupper[ix_['P5']] = v_['PMAX2']
        self.xupper[ix_['P1']] = v_['PMAX1']
        self.xupper[ix_['P7']] = v_['PMAX1']
        self.xupper[ix_['F4']] = 400.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.x0[ix_['P'+str(I)]] = float(965.0)
        self.x0[ix_['Q1']] = float(100.0)
        self.x0[ix_['Q2']] = float(100.0)
        self.x0[ix_['Q3']] = float(-100.0)
        self.x0[ix_['Q4']] = float(-100.0)
        self.x0[ix_['F1']] = float(1000.0)
        self.x0[ix_['F2']] = float(1000.0)
        self.x0[ix_['F3']] = float(1000.0)
        self.x0[ix_['F4']] = float(1000.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eFORQ', iet_)
        elftv = loaset(elftv,it,0,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'PSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            vname = 'P'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'QTO'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eFORQ')
            ielftype = arrset(ielftype, ie, iet_["eFORQ"])
            vname = 'Q'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        v_['K'] = -0.597053452
        ig = ig_['PAN1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['QTO1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['K']))
        ig = ig_['PAN2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['QTO2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['K']))
        ig = ig_['PAN3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ5'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['QTO3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['K']))
        ig = ig_['PAN4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['PSQ7'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['QTO4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['K']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "LOI2-MN-15-11"

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*np.absolute(EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        PO = EV_[0]>0.0
        if PO!=0:
            GO = 2*EV_[0]
        if PO==0:
            GO = -2*EV_[0]
        if PO!=0:
            HO = 2
        if PO==0:
            HO = -2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = GO
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = HO
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eFORQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*np.absolute(EV_[0])**0.8539
        if not isinstance( f_, float ):
            f_   = f_.item();
        POS = EV_[0]>0.0
        if POS!=0:
            GG = 1.8539*EV_[0]**0.8539
        if POS==0:
            GG = 1.8539*np.absolute(EV_[0])**0.8539
        if POS!=0:
            HH = 1.8539*0.8539*EV_[0]**(-0.1461)
        if POS==0:
            HH = -1.8539*0.8539*np.absolute(EV_[0])**(-0.1461)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = GG
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = HH
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

