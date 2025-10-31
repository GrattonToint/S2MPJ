from s2mpjlib import *
class  HANGING(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HANGING
#    *********
# 
#    A catenary problem in 3 dimensions.  A rectangular grid is hung from its
#    4 corners under gravity.  The problem is to determine the resulting shape.
# 
#    Source:  
#    an example in a talk by Nesterova and Vial, LLN, 1994.
# 
#    SIF input: Ph. Toint, November 1994.
# 
#    classification = "C-CLQR2-AY-V-V"
# 
#    dimension of the grid
# 
#           Alternative values for the SIF file parameters:
# IE NX                  3              $-PARAMETER n = 27
# IE NY                  3              $-PARAMETER
# 
# IE NX                  5              $-PARAMETER n = 90
# IE NY                  6              $-PARAMETER
# 
# IE NX                  10             $-PARAMETER n = 300  original value
# IE NY                  10             $-PARAMETER
# 
# IE NX                  20             $-PARAMETER n = 1800
# IE NY                  30             $-PARAMETER
# 
# IE NX                  40             $-PARAMETER n = 3600
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HANGING'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NX'] = int(3);  #  SIF file default value
        else:
            v_['NX'] = int(args[0])
# IE NY                  30             $-PARAMETER
        if nargin<2:
            v_['NY'] = int(3);  #  SIF file default value
        else:
            v_['NY'] = int(args[1])
        v_['LX'] = 1.8
        v_['LY'] = 1.8
        v_['1'] = 1
        v_['NX-1'] = -1+v_['NX']
        v_['NY-1'] = -1+v_['NY']
        v_['LX2'] = v_['LX']*v_['LX']
        v_['LY2'] = v_['LY']*v_['LY']
        v_['RNX'] = float(v_['NX'])
        v_['RNY'] = float(v_['NY'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Y'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Z'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Z'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Z'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                [ig,ig_,_] = s2mpj_ii('RC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'RC'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('DC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'DC'+str(I)+','+str(J))
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
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['RC'+str(I)+','+str(J)],float(v_['LX2'])))
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['DC'+str(I)+','+str(J)],float(v_['LY2'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['Z'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['Z'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['NX']))+','+str(int(v_['1']))]] = v_['RNX']
        self.xupper[ix_['X'+str(int(v_['NX']))+','+str(int(v_['1']))]] = v_['RNX']
        self.xlower[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        self.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['NY']))]] = v_['RNY']
        self.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['NY']))]] = v_['RNY']
        self.xlower[ix_['Z'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        self.xupper[ix_['Z'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNX']
        self.xupper[ix_['X'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNX']
        self.xlower[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNY']
        self.xupper[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNY']
        self.xlower[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = 0.0
        self.xupper[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                self.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['RI-1'])
                self.x0[ix_['Y'+str(I)+','+str(J)]] = float(v_['RJ-1'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for J in range(int(v_['1']),int(v_['NY-1'])+1):
            v_['J+1'] = 1+J
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ename = 'RX'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'RY'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'RZ'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Z'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'DX'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'DY'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'DZ'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eISQ')
                    ielftype = arrset(ielftype,ie,iet_['eISQ'])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Z'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='YY')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                ig = ig_['RC'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RX'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['RY'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RZ'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['DC'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['DX'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['DY'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['DZ'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(3,3)          -6.1184107487
# LO SOLTN(5,6)          -77.260229515
# LO SOLTN(10,10)        -620.17603242
# LO SOLTN(20,30)        -1025.4292887
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQR2-AY-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

