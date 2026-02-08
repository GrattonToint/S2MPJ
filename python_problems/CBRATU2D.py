from s2mpjlib import *
class  CBRATU2D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CBRATU2D
#    *********
# 
#    The complex 2D Bratu problem on the unit square, using finite
#    differences.
# 
#    Source: problem 3 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit square (variable).
#    There are 2*P**2 variables
# 
#           Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 32     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CBRATU2D'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(4);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   7              $-PARAMETER n = 98
# IE P                   16             $-PARAMETER n = 512
# IE P                   23             $-PARAMETER n = 1058
# IE P                   40             $-PARAMETER n = 3200
        if nargin<2:
            v_['LAMBDA'] = float(5.0);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['1.0'] = 1.0
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = v_['1.0']/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['C'] = v_['H2']*v_['LAMBDA']
        v_['-C'] = -1.0*v_['C']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['R'] = 1+I
            v_['S'] = -1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['V'] = 1+J
                v_['W'] = -1+J
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)]])
                valA = np.append(valA,float(4.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(int(v_['R']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(int(v_['S']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(int(v_['V']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(int(v_['W']))]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('F'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'F'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(4.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['R']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['S']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['V']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['W']))]])
                valA = np.append(valA,float(-1.0))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            self.xlower[ix_['U'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xupper[ix_['U'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xlower[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
            self.xupper[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
            self.xlower[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['1']))+','+str(J)]] = 0.0
            self.xlower[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
            self.xupper[ix_['X'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            self.xlower[ix_['U'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xupper[ix_['U'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xlower[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xupper[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xlower[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['P']))]] = 0.0
            self.xlower[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['1']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eRPART', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        [it,iet_,_] = s2mpj_ii( 'eCPART', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eRPART')
                ielftype = arrset(ielftype,ie,iet_["eRPART"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eCPART')
                ielftype = arrset(ielftype,ie,iet_["eCPART"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='V')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-C']))
                ig = ig_['F'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-C']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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
        self.pbclass   = "C-CNOR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eRPART(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPU = np.exp(EV_[0,0])
        EXPUC = EXPU*np.cos(EV_[1,0])
        EXPUS = EXPU*np.sin(EV_[1,0])
        f_   = EXPUC
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPUC
            g_[1] = -EXPUS
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EXPUC
                H_[0,1] = -EXPUS
                H_[1,0] = H_[0,1]
                H_[1,1] = -EXPUC
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCPART(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPU = np.exp(EV_[0,0])
        EXPUC = EXPU*np.cos(EV_[1,0])
        EXPUS = EXPU*np.sin(EV_[1,0])
        f_   = EXPUS
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPUS
            g_[1] = EXPUC
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = EXPUS
                H_[0,1] = EXPUC
                H_[1,0] = H_[0,1]
                H_[1,1] = -EXPUS
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

