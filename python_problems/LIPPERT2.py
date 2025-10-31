from s2mpjlib import *
class  LIPPERT2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LIPPERT2
#    *********
# 
#    A discrete approximation to a continuum optimal flow problem
#    in the unit square. The continuum problem requires that the
#    divergence of a given flow should be given everywhere in the
#    region of interest, with the restriction that the capacity of
#    the flow is bounded. The aim is then to maximize the given flow.
# 
#    The discrete problem (dual formulation 2) in the unit square is to 
#      minimize r
#      subject to dx( u_ij - ui-1j ) + dx( v_ij - vij-1 ) = s_ij
#                 u_ij^2 + v_ij^2 <= r^2
#                 u_i-1j^2 + v_ij^2 <= r^2
#                 u_ij^2 + v_ij-1^2 <= r^2
#                 u_i-1j^2 + v_ij-1^2 <= r^2
#      where 1 <= i <= nx, 1 <= j <= ny
#      and        r >= 0
# 
#    Source: R. A. Lippert
#      "Discrete approximations to continuum optimal flow problems"
#      Tech. Report, Dept of Maths, M.I.T., 2006
#    following a suggestion by Gil Strang
# 
#    SIF input: Nick Gould, September 2006
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CLQR2-MN-V-V"
# 
#    Number of nodes in x direction
# 
#           Alternative values for the SIF file parameters:
# IE NX                  2              $-PARAMETER
# IE NX                  3              $-PARAMETER
# IE NX                  10             $-PARAMETER
# IE NX                  40             $-PARAMETER
# IE NX                  100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LIPPERT2'

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
# IE NY                  2              $-PARAMETER
# IE NY                  3              $-PARAMETER
# IE NY                  10             $-PARAMETER 
# IE NY                  40             $-PARAMETER
# IE NY                  100            $-PARAMETER
        if nargin<2:
            v_['NY'] = int(10);  #  SIF file default value
        else:
            v_['NY'] = int(args[1])
        v_['X+'] = 1+v_['NX']
        v_['X-'] = -1+v_['NX']
        v_['Y+'] = 1+v_['NY']
        v_['Y-'] = -1+v_['NY']
        v_['1'] = 1
        v_['0'] = 0
        v_['HALF'] = 0.5
        v_['ONE'] = 1.0
        v_['-ONE'] = -1.0
        v_['S'] = 1.0
        v_['-S'] = v_['S']*v_['-ONE']
        v_['RX'] = float(v_['NX'])
        v_['DX'] = v_['ONE']/v_['RX']
        v_['-DX'] = v_['-ONE']/v_['RX']
        v_['DX/2'] = v_['DX']*v_['HALF']
        v_['RY'] = float(v_['NY'])
        v_['DY'] = v_['ONE']/v_['RY']
        v_['-DY'] = v_['-ONE']/v_['RY']
        v_['DY/2'] = v_['DY']*v_['HALF']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        self.xnames=arrset(self.xnames,iv,'R')
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'V'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['R']])
        valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'O'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['DX']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['U'+str(int(v_['I-1']))+','+str(J)]])
                valA = np.append(valA,float(v_['-DX']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['DY']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(I)+','+str(int(v_['J-1']))]])
                valA = np.append(valA,float(v_['-DY']))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'B'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'D'+str(I)+','+str(J))
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
            for J in range(int(v_['1']),int(v_['NY'])+1):
                self.gconst = arrset(self.gconst,ig_['O'+str(I)+','+str(J)],float(v_['S']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['R']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['R']] = float(1.0)
        v_['ALPHA'] = 0.0
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                self.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['ALPHA'])
            v_['ALPHA'] = v_['ALPHA']+v_['DX/2']
        v_['ALPHA'] = 0.0
        for J in range(int(v_['0']),int(v_['NY'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                self.x0[ix_['V'+str(I)+','+str(J)]] = float(v_['ALPHA'])
            v_['ALPHA'] = v_['ALPHA']+v_['DX/2']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'ALPHA')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'RHO2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='ALPHA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eSQR')
                    ielftype = arrset(ielftype,ie,iet_['eSQR'])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='ALPHA')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                ename = 'Q'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    self.elftype = arrset(self.elftype,ie,'eSQR')
                    ielftype = arrset(ielftype,ie,iet_['eSQR'])
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='ALPHA')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                ig = ig_['A'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['Q'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['B'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['P'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['Q'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['C'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['Q'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['D'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['P'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['Q'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               3.77245385
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

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
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

