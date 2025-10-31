from s2mpjlib import *
class  OPTMASS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTMASS
#    *********
# 
#    A constrained optimal control problem
#    adapted from Gawande and Dunn
# 
#    The problem is that of a particle of unit mass moving on a
#    frictionless plane under the action of a controlling force whose
#    magnitude may not exceed unity. At time=0, the particle moves through
#    the origin of the plane in the direction of the positive x-axis with
#    speed SPEED.  The cost function incorporates two conflicting control
#    objectives, namely: maximization of the particle's final (at time=1)
#    distance from the origin and minimization of its final speed.  By
#    increasing the  value of the penalty constant PEN, more stress can be
#    placed on the latter objective.
# 
#    Gawande and Dunn originally use a starting point (in the control
#    only) that is much closer to the solution than the one chosen
#    here.
# 
#    Source:
#    M. Gawande and J. Dunn,
#    "A Projected Newton Method in a Cartesian Product of Balls",
#    JOTA 59(1): 59-69, 1988.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-CQQR2-AN-V-V"
# 
#    Number of discretization steps in the time interval
#    The number of variables is 6 * (N + 2) -2 , 4 of which are fixed.
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n = 70    original value
# IE N                   100            $-PARAMETER n = 610
# IE N                   200            $-PARAMETER n = 1210
# IE N                   500            $-PARAMETER n = 3010
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OPTMASS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   1000           $-PARAMETER n = 6010
# IE N                   5000           $-PARAMETER n = 30010
        v_['SPEED'] = 0.01
        v_['PEN'] = 0.335
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N+1'] = 1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['1/N'] = 1.0/v_['RN']
        v_['-1/N'] = -1.0*v_['1/N']
        v_['1/N2'] = v_['1/N']*v_['1/N']
        v_['-1/2N2'] = -0.5*v_['1/N2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['2'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(J)+','+str(I))
                [iv,ix_,_] = s2mpj_ii('V'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'V'+str(J)+','+str(I))
                [iv,ix_,_] = s2mpj_ii('F'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'F'+str(J)+','+str(I))
        for J in range(int(v_['1']),int(v_['2'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(int(v_['N+1'])),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(J)+','+str(int(v_['N+1'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(J)+','+str(int(v_['N+1'])),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(J)+','+str(int(v_['N+1'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('F',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N+1'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['2'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'A'+str(J)+','+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)+','+str(I)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(J)+','+str(int(v_['I-1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(J)+','+str(int(v_['I-1']))]])
                valA = np.append(valA,float(v_['-1/N']))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['F'+str(J)+','+str(int(v_['I-1']))]])
                valA = np.append(valA,float(v_['-1/2N2']))
                [ig,ig_,_] = s2mpj_ii('B'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B'+str(J)+','+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(J)+','+str(I)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['V'+str(J)+','+str(int(v_['I-1']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['F'+str(J)+','+str(int(v_['I-1']))]])
                valA = np.append(valA,float(v_['-1/N']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(I))
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
        for I in range(int(v_['0']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['C'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]] = v_['SPEED']
        self.xupper[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]] = v_['SPEED']
        self.xlower[ix_['V'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['V'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        if('V'+str(int(v_['1']))+','+str(int(v_['0'])) in ix_):
            self.x0[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]]  = (
                  float(v_['SPEED']))
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['V'+str(int(v_['1']))+','+str(int(v_['0']))]),float(v_['SPEED'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'O1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['N+1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['N+1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'V'+str(int(v_['1']))+','+str(int(v_['N+1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'V'+str(int(v_['2']))+','+str(int(v_['N+1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['2'])+1):
                ename = 'D'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype,ie,iet_["eSQ"])
                vname = 'F'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['F']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['O2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['PEN']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['PEN']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['C'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['1']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['2']))+','+str(I)]))
            self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           -0.04647
# LO SOLTN(100)          ???
# LO SOLTN(200)          ???
# LO SOLTN(500)          ???
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQQR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

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

