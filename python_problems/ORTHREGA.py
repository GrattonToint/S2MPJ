from s2mpjlib import *
class  ORTHREGA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGA
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) an ellipse to a set of points
#    in the plane.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-CQQR2-AN-V-V"
# 
#    Number of levels in the generation of the data points
#    ( number of data points =     4**LEVELS
#      number of variables   = 2 * 4**LEVELS + 5
#      number of constraints =     4**LEVELS         )
# 
#           Alternative values for the SIF file parameters:
# IE LEVELS              3              $-PARAMETER n = 133    original value
# IE LEVELS              4              $-PARAMETER n = 517
# IE LEVELS              5              $-PARAMETER n = 2053
# IE LEVELS              6              $-PARAMETER n = 8197
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHREGA'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['LEVELS'] = int(2);  #  SIF file default value
        else:
            v_['LEVELS'] = int(args[0])
# IE LEVELS              7              $-PARAMETER n = 32773
# IE LEVELS              8              $-PARAMETER n = 131077
        v_['A'] = 9.0
        v_['B'] = 6.0
        v_['CX'] = 0.5
        v_['CY'] = 0.5
        v_['1'] = 1
        v_['PI'] = 3.1415926535
        v_['-A'] = -1.0*v_['A']
        v_['-B'] = -1.0*v_['B']
        v_['NPTS'] = 1
        v_['XD'+str(int(v_['1']))] = v_['CX']
        v_['YD'+str(int(v_['1']))] = v_['CY']
        for I in range(int(v_['1']),int(v_['LEVELS'])+1):
            v_['NP'] = 0+v_['NPTS']
            for J in range(int(v_['1']),int(v_['NP'])+1):
                v_['XZ'+str(J)] = v_['XD'+str(J)]
                v_['YZ'+str(J)] = v_['YD'+str(J)]
            v_['NPTS'] = 0
            for J in range(int(v_['1']),int(v_['NP'])+1):
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['A']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['A']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['B']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['-B']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['-A']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['-A']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['-B']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['B']
            v_['A'] = v_['A']/v_['PI']
            v_['B'] = v_['B']/v_['PI']
            v_['-A'] = v_['-A']/v_['PI']
            v_['-B'] = v_['-B']/v_['PI']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('H11',ix_)
        self.xnames=arrset(self.xnames,iv,'H11')
        [iv,ix_,_] = s2mpj_ii('H12',ix_)
        self.xnames=arrset(self.xnames,iv,'H12')
        [iv,ix_,_] = s2mpj_ii('H22',ix_)
        self.xnames=arrset(self.xnames,iv,'H22')
        [iv,ix_,_] = s2mpj_ii('G1',ix_)
        self.xnames=arrset(self.xnames,iv,'G1')
        [iv,ix_,_] = s2mpj_ii('G2',ix_)
        self.xnames=arrset(self.xnames,iv,'G2')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [ig,ig_,_] = s2mpj_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('OY'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Y'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            self.gconst = arrset(self.gconst,ig_['OX'+str(I)],float(v_['XD'+str(I)]))
            self.gconst = arrset(self.gconst,ig_['OY'+str(I)],float(v_['YD'+str(I)]))
            self.gconst = arrset(self.gconst,ig_['E'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('H11' in ix_):
            self.x0[ix_['H11']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['H11']),float(1.0)))
        if('H12' in ix_):
            self.x0[ix_['H12']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['H12']),float(0.0)))
        if('H22' in ix_):
            self.x0[ix_['H22']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['H22']),float(1.0)))
        if('G1' in ix_):
            self.x0[ix_['G1']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['G1']),float(0.0)))
        if('G2' in ix_):
            self.x0[ix_['G2']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['G2']),float(0.0)))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            if('X'+str(I) in ix_):
                self.x0[ix_['X'+str(I)]] = float(v_['XD'+str(I)])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['XD'+str(I)])))
            if('Y'+str(I) in ix_):
                self.x0[ix_['Y'+str(I)]] = float(v_['YD'+str(I)])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y'+str(I)]),float(v_['YD'+str(I)])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eHXX', iet_)
        elftv = loaset(elftv,it,0,'H')
        elftv = loaset(elftv,it,1,'X')
        [it,iet_,_] = s2mpj_ii( 'eHXY', iet_)
        elftv = loaset(elftv,it,0,'H')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        [it,iet_,_] = s2mpj_ii( 'eGX', iet_)
        elftv = loaset(elftv,it,0,'G')
        elftv = loaset(elftv,it,1,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eHXX')
            ielftype = arrset(ielftype,ie,iet_["eHXX"])
            vname = 'H11'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='H')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eHXY')
            ielftype = arrset(ielftype,ie,iet_["eHXY"])
            vname = 'H12'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='H')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eHXX')
            ielftype = arrset(ielftype,ie,iet_["eHXX"])
            vname = 'H22'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='H')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'ED'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eGX')
            ielftype = arrset(ielftype,ie,iet_["eGX"])
            vname = 'G1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='G')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EE'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eGX')
            ielftype = arrset(ielftype,ie,iet_["eGX"])
            vname = 'G2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='G')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ig = ig_['OX'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['OY'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(2.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['ED'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-2.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
# LO SOLTN(3)             350.29936756
# LO SOLTN(4)             1414.0524915
# LO SOLTN(5)             ???
# LO SOLTN(6)             ???
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
        self.pbclass   = "C-CQQR2-AN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHXX(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[1,0]
            g_[1] = 2.0*EV_[0,0]*EV_[1,0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EV_[1,0]+EV_[1,0]
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0,0]+EV_[0,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eHXY(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]
            g_[1] = EV_[0,0]*EV_[2,0]
            g_[2] = EV_[0,0]*EV_[1,0]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0,0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eGX(self, nargout,*args):

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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
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

