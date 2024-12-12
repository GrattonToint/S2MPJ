from s2mpjlib import *
class  PENTAGON(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PENTAGON
#    *********
# 
#    An approximation to the problem of finding 3 points in a 2D
#    pentagon whose minimal distance is maximal.
# 
#    Source:
#    M.J.D. Powell,
#    " TOLMIN: a Fortran package for linearly constrained
#    optimization problems",
#    Report DAMTP 1989/NA2, University of Cambridge, UK, 1989.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "C-COLR2-AY-6-15"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PENTAGON'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['3'] = 3
        v_['4'] = 4
        v_['2PI/5'] = 1.2566371
        for J in range(int(v_['0']),int(v_['4'])+1):
            v_['RJ'] = float(J)
            v_['TJ'] = v_['2PI/5']*v_['RJ']
            v_['C'+str(J)] = np.cos(v_['TJ'])
            v_['S'+str(J)] = np.sin(v_['TJ'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['3'])+1):
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
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['3'])+1):
            for J in range(int(v_['0']),int(v_['4'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)]])
                valA = np.append(valA,float(v_['C'+str(J)]))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)]])
                valA = np.append(valA,float(v_['S'+str(J)]))
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),1.0)
        self.gconst = arrset(self.gconst,ig_['OBJ'],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(-1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(-1.0)))
        if('Y1' in ix_):
            self.x0[ix_['Y1']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y1']),float(0.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(0.0)))
        if('Y2' in ix_):
            self.x0[ix_['Y2']] = float(-1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y2']),float(-1.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(1.0)))
        if('Y3' in ix_):
            self.x0[ix_['Y3']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['Y3']),float(1.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eIDIST', iet_)
        elftv = loaset(elftv,it,0,'XA')
        elftv = loaset(elftv,it,1,'YA')
        elftv = loaset(elftv,it,2,'XB')
        elftv = loaset(elftv,it,3,'YB')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'D12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype,ie,iet_["eIDIST"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype,ie,iet_["eIDIST"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eIDIST')
        ielftype = arrset(ielftype,ie,iet_["eIDIST"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YA')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='XB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'Y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='YB')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['D13'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D32'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
#    Solution
#  LO SOLTN              1.36521631D-04
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
        self.pbclass   = "C-COLR2-AY-6-15"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eIDIST(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        D = IV_[0]*IV_[0]+IV_[1]*IV_[1]
        D9 = D**9
        D10 = D9*D
        f_   = 1.0/D**8
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -16.0*IV_[0]/D9
            g_[1] = -16.0*IV_[1]/D9
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 16.0*(18.0*IV_[0]*IV_[0]-D)/D10
                H_[0,1] = 288.0*IV_[0]*IV_[1]/D10
                H_[1,0] = H_[0,1]
                H_[1,1] = 16.0*(18.0*IV_[1]*IV_[1]-D)/D10
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

