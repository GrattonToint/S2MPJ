from s2mpjlib import *
class  HS70(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS70
#    *********
# 
#    This problem arises in water flow routing.
# 
#    Source: problem 70 incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991, modified May 2024
# 
#    classification = "C-CSQR2-MN-4-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS70'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
        v_['1'] = 1
        v_['2'] = 2
        v_['19'] = 19
        v_['C1'] = 0.1
        v_['C2'] = 1.0
        v_['C3'] = 2.0
        v_['C4'] = 3.0
        v_['C5'] = 4.0
        v_['C6'] = 5.0
        v_['C7'] = 6.0
        v_['C8'] = 7.0
        v_['C9'] = 8.0
        v_['C10'] = 9.0
        v_['C11'] = 10.0
        v_['C12'] = 11.0
        v_['C13'] = 12.0
        v_['C14'] = 13.0
        v_['C15'] = 14.0
        v_['C16'] = 15.0
        v_['C17'] = 16.0
        v_['C18'] = 17.0
        v_['C19'] = 18.0
        v_['Y1'] = 0.00189
        v_['Y2'] = 0.1038
        v_['Y3'] = 0.268
        v_['Y4'] = 0.506
        v_['Y5'] = 0.577
        v_['Y6'] = 0.604
        v_['Y7'] = 0.725
        v_['Y8'] = 0.898
        v_['Y9'] = 0.947
        v_['Y10'] = 0.845
        v_['Y11'] = 0.702
        v_['Y12'] = 0.528
        v_['Y13'] = 0.385
        v_['Y14'] = 0.257
        v_['Y15'] = 0.159
        v_['Y16'] = 0.0869
        v_['Y17'] = 0.0453
        v_['Y18'] = 0.01509
        v_['Y19'] = 0.00189
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['19'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0e+0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0e+0))
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
        for I in range(int(v_['1']),int(v_['19'])+1):
            self.gconst = arrset(self.gconst,ig_['OBJ'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.00001)
        self.xupper = np.full((self.n,1),100.0)
        self.xupper[ix_['X3']] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(2.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(2.0)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(4.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(4.0))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(0.04)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(0.04)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(2.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(2.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eY1', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'eY2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X3')
        elftv = loaset(elftv,it,1,'X4')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['19'])+1):
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eY1')
            ielftype = arrset(ielftype,ie,iet_["eY1"])
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='C')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(I)]))
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eY2')
            ielftype = arrset(ielftype,ie,iet_["eY2"])
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='C')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(I)]))
        ename = 'C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
        posev = np.where(elftv[ielftype[ie]]=='X3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.00001),float(100.0),None)
        posev = np.where(elftv[ielftype[ie]]=='X4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQR',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['19'])+1):
            ig = ig_['OBJ'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gSQR')
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['Y'+str(I)+','+str(int(v_['1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['Y'+str(I)+','+str(int(v_['2']))]))
            self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.007498464
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSQR2-MN-4-1"
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
        self.efpar = arrset( self.efpar,0,np.sqrt(1.0e+0/6.2832e+0))
        return pbm

    @staticmethod
    def eY1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B = EV_[1]+EV_[2]*(1.0e+0-EV_[1])
        CI = self.elpar[iel_][0]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[0])
        P0V1 = -1.0e+0/(1.2e+1*EV_[0]**2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[0]**3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0**2)
        P1V1V1 = (2.0e+0*P0V1**2/P0-P0V1V1)/(P0**2)
        P2 = EV_[1]
        P3 = B**EV_[0]
        LOGB = np.log(B)
        P3V1 = P3*LOGB
        P3V2 = EV_[0]*(1.0e+0-EV_[2])*B**(EV_[0]-1.0e+0)
        P3V3 = EV_[0]*(1.0e+0-EV_[1])*B**(EV_[0]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*(1.0e+0-EV_[2])/B
        P3V1V3 = P3V3*LOGB+P3*(1.0e+0-EV_[1])/B
        P3V2V2 = EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[2])**2*B**(EV_[0]-1.0e+0)
        P3V2V3  = (
              -EV_[0]*B**(EV_[0]-1.0e+0)+EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[1])*(1.0e+0-EV_[2])*B**(EV_[0]-2.0e+0))
        P3V3V3 = EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[1])**2*B**(EV_[0]-2.0e+0)
        P4 = self.efpar[0]*np.sqrt(EV_[0])
        P4V1 = 5.0e-1*self.efpar[0]*np.sqrt(1.0e+0/EV_[0])
        P4V1V1 = -2.5e-1*self.efpar[0]*np.sqrt(1.0e+0/EV_[0]**3)
        C5 = CI**(-1.0e+0)
        P5 = C5*CI**EV_[0]
        P5V1 = P5*np.log(CI)
        P5V1V1 = P5V1*np.log(CI)
        P6 = np.exp(EV_[0]*(1.0e+0-CI*B))
        P6V1 = P6*(1.0e+0-CI*B)
        P6V2 = -P6*EV_[0]*CI*(1.0e+0-EV_[2])
        P6V3 = -P6*EV_[0]*CI*(1.0e+0-EV_[1])
        P6V1V1 = P6*(1.0e+0-CI*B)**2
        P6V1V2 = P6V2*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[2])
        P6V1V3 = P6V3*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[1])
        P6V2V2 = -P6V2*EV_[0]*CI*(1.0e+0-EV_[2])
        P6V2V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[2])+P6*EV_[0]*CI
        P6V3V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[1])
        f_   = P1*P2*P3*P4*P5*P6
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[1] = P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[2] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[0,1]  = (
                      P1V1*(P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(P3V1*P4*P5*P6+(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[2,0] = H_[0,2]
                H_[1,1]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2+P3*P6V2+P3V2*P6)))
                H_[1,2]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3+P3V3*P6+P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[2,1] = H_[1,2]
                H_[2,2] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B = EV_[1]+EV_[2]*(1.0e+0-EV_[1])
        CI = self.elpar[iel_][0]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[0])
        P0V1 = -1.0e+0/(1.2e+1*EV_[0]**2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[0]**3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0**2)
        P1V1V1 = (2.0e+0*P0V1**2/P0-P0V1V1)/(P0**2)
        P2 = 1.0e+0-EV_[1]
        P3 = (B/EV_[2])**EV_[0]
        LOGB = np.log(B/EV_[2])
        P3V1 = P3*LOGB
        P3V2 = EV_[0]*(-1.0e+0+1.0e+0/EV_[2])*(B/EV_[2])**(EV_[0]-1.0e+0)
        P3V3 = -EV_[0]*(EV_[1]/EV_[2]**2)*(B/EV_[2])**(EV_[0]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*EV_[2]*(-1.0e+0+1.0e+0/EV_[2])/B
        P3V1V3 = P3V3*LOGB-P3*EV_[1]/(B*EV_[2])
        P3V2V2  = (
              EV_[0]*(EV_[0]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[2])**2*(B/EV_[2])**(EV_[0]-2.0e+0))
        P3V2V3  = (
              EV_[0]*(-1.0e+0/EV_[2]**2)*(B/EV_[2])**(EV_[0]-1.0e+0)+EV_[0]*(EV_[0]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[2])*(-EV_[1]/EV_[2]**2)*(B/EV_[2])**(EV_[0]-2.0e+0))
        P3V3V3 = (2.0e+0*EV_[0]*(EV_[1]/EV_[2]**3)*(B/EV_[2])**(EV_[0]-1.0e+0)+
             EV_[0]*(EV_[0]-1.0e+0)*(EV_[1]/EV_[2]**2)**2*(B/EV_[2])**(EV_[0]-2.0e+0))
        P4 = self.efpar[0]*np.sqrt(EV_[0])
        P4V1 = 5.0e-1*self.efpar[0]*np.sqrt(1.0e+0/EV_[0])
        P4V1V1 = -2.5e-1*self.efpar[0]*np.sqrt(1.0e+0/EV_[0]**3)
        C5 = CI**(-1.0e+0)
        P5 = C5*CI**EV_[0]
        P5V1 = P5*np.log(CI)
        P5V1V1 = P5V1*np.log(CI)
        P6 = np.exp(EV_[0]*(1.0e+0-CI*B/EV_[2]))
        P6V1 = P6*(1.0e+0-CI*B/EV_[2])
        P6V2 = -P6*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]
        P6V3 = P6*EV_[0]*CI*EV_[1]/EV_[2]**2
        P6V1V1 = P6*(1.0e+0-CI*B/EV_[2])**2
        P6V1V2 = P6V2*(1.0e+0-CI*B/EV_[2])-P6*CI*(-1.0e+0+1.0e+0/EV_[2])
        P6V1V3 = P6V3*(1.0e+0-CI*B/EV_[2])+P6*CI*EV_[1]/EV_[2]**2
        P6V2V2 = -P6V2*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]
        P6V2V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]+P6*EV_[0]*CI/EV_[2]**2
        P6V3V3 = (P6V3*EV_[0]*CI*EV_[1]/EV_[2]**2-2.0e+0*P6*EV_[0]*CI*EV_[1]/
             EV_[2]**3)
        f_   = P1*P2*P3*P4*P5*P6
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[1] = -P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[2] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[0,1]  = (
                      P1V1*(-P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(-P3V1*P4*P5*P6-(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[2,0] = H_[0,2]
                H_[1,1]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2-P3*P6V2-P3V2*P6)))
                H_[1,2]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3-P3V3*P6-P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[2,1] = H_[1,2]
                H_[2,2] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
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
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0e+0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQR(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

