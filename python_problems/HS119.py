from s2mpjlib import *
class  HS119(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS119
#    *********
# 
#    Source: problem 119 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 7 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-COLR2-AN-16-8"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS119'

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
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['15'] = 15
        v_['16'] = 16
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
        for I in range(int(v_['1']),int(v_['8'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                v_['B'+str(I)+','+str(J)] = 0.0
        for I in range(int(v_['1']),int(v_['16'])+1):
            v_['A'+str(I)+','+str(I)] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['4']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['8']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['3']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['9']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['11']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['6']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['12']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['6']))+','+str(int(v_['8']))] = 1.0
        v_['A'+str(int(v_['6']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['7']))+','+str(int(v_['11']))] = 1.0
        v_['A'+str(int(v_['7']))+','+str(int(v_['13']))] = 1.0
        v_['A'+str(int(v_['8']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['8']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['9']))+','+str(int(v_['12']))] = 1.0
        v_['A'+str(int(v_['9']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['10']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['11']))+','+str(int(v_['13']))] = 1.0
        v_['A'+str(int(v_['12']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['13']))+','+str(int(v_['14']))] = 1.0
        v_['B'+str(int(v_['1']))+','+str(int(v_['1']))] = 0.22
        v_['B'+str(int(v_['2']))+','+str(int(v_['1']))] = -1.46
        v_['B'+str(int(v_['3']))+','+str(int(v_['1']))] = 1.29
        v_['B'+str(int(v_['4']))+','+str(int(v_['1']))] = -1.10
        v_['B'+str(int(v_['7']))+','+str(int(v_['1']))] = 1.12
        v_['B'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.20
        v_['B'+str(int(v_['3']))+','+str(int(v_['2']))] = -0.89
        v_['B'+str(int(v_['4']))+','+str(int(v_['2']))] = -1.06
        v_['B'+str(int(v_['6']))+','+str(int(v_['2']))] = -1.72
        v_['B'+str(int(v_['8']))+','+str(int(v_['2']))] = 0.45
        v_['B'+str(int(v_['1']))+','+str(int(v_['3']))] = 0.19
        v_['B'+str(int(v_['2']))+','+str(int(v_['3']))] = -1.30
        v_['B'+str(int(v_['4']))+','+str(int(v_['3']))] = 0.95
        v_['B'+str(int(v_['6']))+','+str(int(v_['3']))] = -0.33
        v_['B'+str(int(v_['8']))+','+str(int(v_['3']))] = 0.26
        v_['B'+str(int(v_['1']))+','+str(int(v_['4']))] = 0.25
        v_['B'+str(int(v_['2']))+','+str(int(v_['4']))] = 1.82
        v_['B'+str(int(v_['4']))+','+str(int(v_['4']))] = -0.54
        v_['B'+str(int(v_['5']))+','+str(int(v_['4']))] = -1.43
        v_['B'+str(int(v_['7']))+','+str(int(v_['4']))] = 0.31
        v_['B'+str(int(v_['8']))+','+str(int(v_['4']))] = -1.10
        v_['B'+str(int(v_['1']))+','+str(int(v_['5']))] = 0.15
        v_['B'+str(int(v_['2']))+','+str(int(v_['5']))] = -1.15
        v_['B'+str(int(v_['3']))+','+str(int(v_['5']))] = -1.16
        v_['B'+str(int(v_['5']))+','+str(int(v_['5']))] = 1.51
        v_['B'+str(int(v_['6']))+','+str(int(v_['5']))] = 1.62
        v_['B'+str(int(v_['8']))+','+str(int(v_['5']))] = 0.58
        v_['B'+str(int(v_['1']))+','+str(int(v_['6']))] = 0.11
        v_['B'+str(int(v_['3']))+','+str(int(v_['6']))] = -0.96
        v_['B'+str(int(v_['4']))+','+str(int(v_['6']))] = -1.78
        v_['B'+str(int(v_['5']))+','+str(int(v_['6']))] = 0.59
        v_['B'+str(int(v_['6']))+','+str(int(v_['6']))] = 1.24
        v_['B'+str(int(v_['1']))+','+str(int(v_['7']))] = 0.12
        v_['B'+str(int(v_['2']))+','+str(int(v_['7']))] = 0.80
        v_['B'+str(int(v_['4']))+','+str(int(v_['7']))] = -0.41
        v_['B'+str(int(v_['5']))+','+str(int(v_['7']))] = -0.33
        v_['B'+str(int(v_['6']))+','+str(int(v_['7']))] = 0.21
        v_['B'+str(int(v_['7']))+','+str(int(v_['7']))] = 1.12
        v_['B'+str(int(v_['8']))+','+str(int(v_['7']))] = -1.03
        v_['B'+str(int(v_['1']))+','+str(int(v_['8']))] = 0.13
        v_['B'+str(int(v_['3']))+','+str(int(v_['8']))] = -0.49
        v_['B'+str(int(v_['5']))+','+str(int(v_['8']))] = -0.43
        v_['B'+str(int(v_['6']))+','+str(int(v_['8']))] = -0.26
        v_['B'+str(int(v_['8']))+','+str(int(v_['8']))] = 0.10
        v_['B'+str(int(v_['1']))+','+str(int(v_['9']))] = 1.00
        v_['B'+str(int(v_['7']))+','+str(int(v_['9']))] = -0.36
        v_['B'+str(int(v_['2']))+','+str(int(v_['10']))] = 1.00
        v_['B'+str(int(v_['3']))+','+str(int(v_['11']))] = 1.00
        v_['B'+str(int(v_['4']))+','+str(int(v_['12']))] = 1.00
        v_['B'+str(int(v_['5']))+','+str(int(v_['13']))] = 1.00
        v_['B'+str(int(v_['6']))+','+str(int(v_['14']))] = 1.00
        v_['B'+str(int(v_['7']))+','+str(int(v_['15']))] = 1.00
        v_['B'+str(int(v_['8']))+','+str(int(v_['16']))] = 1.00
        v_['C'+str(int(v_['1']))] = 2.5
        v_['C'+str(int(v_['2']))] = 1.1
        v_['C'+str(int(v_['3']))] = -3.1
        v_['C'+str(int(v_['4']))] = -3.5
        v_['C'+str(int(v_['5']))] = 1.3
        v_['C'+str(int(v_['6']))] = 2.1
        v_['C'+str(int(v_['7']))] = 2.3
        v_['C'+str(int(v_['8']))] = -1.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['16'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            [ig,ig_,_] = s2mpj_ii('OG'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['8'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)]])
                valA = np.append(valA,float(v_['B'+str(I)+','+str(J)]))
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
        for I in range(int(v_['1']),int(v_['8'])+1):
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),5.0)
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(10.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftp = []
        elftp = loaset(elftp,it,0,'AIJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'ePROD')
                    ielftype = arrset(ielftype,ie,iet_['ePROD'])
                vname = 'X'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(5.0),float(10.0))
                posev = np.where(elftv[ielftype[ie]]=='U1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(5.0),float(10.0))
                posev = np.where(elftv[ielftype[ie]]=='U2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='AIJ')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['A'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                ig = ig_['OG'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
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
        self.pbclass   = "C-COLR2-AN-16-8"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TU1P1 = 2.0*EV_[0]+1
        TU2P1 = 2.0*EV_[1]+1
        FIRST = EV_[0]**2+EV_[0]+1.0
        SECOND = EV_[1]**2+EV_[1]+1.0
        f_   = self.elpar[iel_][0]*FIRST*SECOND
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*TU1P1*SECOND
            g_[1] = self.elpar[iel_][0]*TU2P1*FIRST
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = self.elpar[iel_][0]*2.0*SECOND
                H_[0,1] = self.elpar[iel_][0]*TU1P1*TU2P1
                H_[1,0] = H_[0,1]
                H_[1,1] = self.elpar[iel_][0]*2.0*FIRST
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

