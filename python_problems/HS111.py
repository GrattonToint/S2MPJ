from s2mpjlib import *
class  HS111(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS111
#    *********
# 
#    This problem is a chemical equilibrium problem involving 3 linear
#    equality constraints.
# 
#    Source: problem 111 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-COOR2-AN-10-3"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS111'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 10
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
        v_['C1'] = -6.089
        v_['C2'] = -17.164
        v_['C3'] = -34.054
        v_['C4'] = -5.914
        v_['C5'] = -24.721
        v_['C6'] = -14.986
        v_['C7'] = -24.100
        v_['C8'] = -10.708
        v_['C9'] = -26.662
        v_['C10'] = -22.179
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
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON2')
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON3')
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
        self.gconst = arrset(self.gconst,ig_['CON1'],float(2.0))
        self.gconst = arrset(self.gconst,ig_['CON2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CON3'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-100.0)
        self.xupper = np.full((self.n,1),100.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(-2.3))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBJ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        elftp = []
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'+str(I)] = float(I)
            v_['RI'+str(I)] = 0.1+v_['RI'+str(I)]
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'O'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eOBJ')
            ielftype = arrset(ielftype,ie,iet_["eOBJ"])
            v_['TEMP'] = v_['RI'+str(int(v_['1']))]
            v_['RI'+str(int(v_['1']))] = v_['RI'+str(I)]
            v_['RI'+str(I)] = v_['TEMP']
            v_['R'] = v_['RI'+str(int(v_['1']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['2']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['3']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['4']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['5']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['6']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['7']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V7')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['8']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['9']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V9')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['10']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='V10')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='C')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['C'+str(I)]))
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eEXP')
            ielftype = arrset(ielftype,ie,iet_["eEXP"])
            vname = 'X'+str(I)
            [iv,ix_]  = (
                  s2mpj_nlx(self,vname,ix_,1,float(-100.0),float(100.0),float(-2.3)))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['O'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CON1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(2.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['CON2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['CON3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -47.707579
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AN-10-3"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EX = np.exp(EV_[0,0])
        f_   = EX
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EX
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOBJ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E1 = np.exp(EV_[0,0])
        E2 = np.exp(EV_[1,0])
        E3 = np.exp(EV_[2,0])
        E4 = np.exp(EV_[3,0])
        E5 = np.exp(EV_[4,0])
        E6 = np.exp(EV_[5,0])
        E7 = np.exp(EV_[6,0])
        E8 = np.exp(EV_[7,0])
        E9 = np.exp(EV_[8,0])
        E10 = np.exp(EV_[9,0])
        SUM = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10
        f_   = E1*(self.elpar[iel_][0]+EV_[0,0]-np.log(SUM))
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E1*(self.elpar[iel_][0]+EV_[0,0]-np.log(SUM))+E1*(1.0e+0-E1/SUM)
            g_[1] = -E1*E2/SUM
            g_[2] = -E1*E3/SUM
            g_[3] = -E1*E4/SUM
            g_[4] = -E1*E5/SUM
            g_[5] = -E1*E6/SUM
            g_[6] = -E1*E7/SUM
            g_[7] = -E1*E8/SUM
            g_[8] = -E1*E9/SUM
            g_[9] = -E1*E10/SUM
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,0] = (E1*(self.elpar[iel_][0]+EV_[0,0]-np.log(SUM))+E1*(1.0e+0-E1/SUM)+
                     E1*(1.0e+0-E1/SUM)+E1*(-E1/SUM)+E1*(E1**2/SUM**2))
                H_[0,1] = (-1.0e+0+E1/SUM)*E1*E2/SUM
                H_[1,0] = H_[0,1]
                H_[1,1] = (-1.0e+0+E2/SUM)*E1*E2/SUM
                H_[0,2] = (-1.0e+0+E1/SUM)*E1*E3/SUM
                H_[2,0] = H_[0,2]
                H_[1,2] = E1*E2*E3/SUM**2
                H_[2,1] = H_[1,2]
                H_[2,2] = (-1.0e+0+E3/SUM)*E1*E3/SUM
                H_[0,3] = (-1.0e+0+E1/SUM)*E1*E4/SUM
                H_[3,0] = H_[0,3]
                H_[1,3] = E1*E2*E4/SUM**2
                H_[3,1] = H_[1,3]
                H_[2,3] = E1*E3*E4/SUM**2
                H_[3,2] = H_[2,3]
                H_[3,3] = (-1.0e+0+E4/SUM)*E1*E4/SUM
                H_[0,4] = (-1.0e+0+E1/SUM)*E1*E5/SUM
                H_[4,0] = H_[0,4]
                H_[1,4] = E1*E2*E5/SUM**2
                H_[4,1] = H_[1,4]
                H_[2,4] = E1*E3*E5/SUM**2
                H_[4,2] = H_[2,4]
                H_[3,4] = E1*E4*E5/SUM**2
                H_[4,3] = H_[3,4]
                H_[4,4] = (-1.0e+0+E5/SUM)*E1*E5/SUM
                H_[0,5] = (-1.0e+0+E1/SUM)*E1*E6/SUM
                H_[5,0] = H_[0,5]
                H_[1,5] = E1*E2*E6/SUM**2
                H_[5,1] = H_[1,5]
                H_[2,5] = E1*E3*E6/SUM**2
                H_[5,2] = H_[2,5]
                H_[3,5] = E1*E4*E6/SUM**2
                H_[5,3] = H_[3,5]
                H_[4,5] = E1*E5*E6/SUM**2
                H_[5,4] = H_[4,5]
                H_[5,5] = (-1.0e+0+E6/SUM)*E1*E6/SUM
                H_[0,6] = (-1.0e+0+E1/SUM)*E1*E7/SUM
                H_[6,0] = H_[0,6]
                H_[1,6] = E1*E2*E7/SUM**2
                H_[6,1] = H_[1,6]
                H_[2,6] = E1*E3*E7/SUM**2
                H_[6,2] = H_[2,6]
                H_[3,6] = E1*E4*E7/SUM**2
                H_[6,3] = H_[3,6]
                H_[4,6] = E1*E5*E7/SUM**2
                H_[6,4] = H_[4,6]
                H_[5,6] = E1*E6*E7/SUM**2
                H_[6,5] = H_[5,6]
                H_[6,6] = (-1.0e+0+E7/SUM)*E1*E7/SUM
                H_[0,7] = (-1.0e+0+E1/SUM)*E1*E8/SUM
                H_[7,0] = H_[0,7]
                H_[1,7] = E1*E2*E8/SUM**2
                H_[7,1] = H_[1,7]
                H_[2,7] = E1*E3*E8/SUM**2
                H_[7,2] = H_[2,7]
                H_[3,7] = E1*E4*E8/SUM**2
                H_[7,3] = H_[3,7]
                H_[4,7] = E1*E5*E8/SUM**2
                H_[7,4] = H_[4,7]
                H_[5,7] = E1*E6*E8/SUM**2
                H_[7,5] = H_[5,7]
                H_[6,7] = E1*E7*E8/SUM**2
                H_[7,6] = H_[6,7]
                H_[7,7] = (-1.0e+0+E8/SUM)*E1*E8/SUM
                H_[0,8] = (-1.0e+0+E1/SUM)*E1*E9/SUM
                H_[8,0] = H_[0,8]
                H_[1,8] = E1*E2*E9/SUM**2
                H_[8,1] = H_[1,8]
                H_[2,8] = E1*E3*E9/SUM**2
                H_[8,2] = H_[2,8]
                H_[3,8] = E1*E4*E9/SUM**2
                H_[8,3] = H_[3,8]
                H_[4,8] = E1*E5*E9/SUM**2
                H_[8,4] = H_[4,8]
                H_[5,8] = E1*E6*E9/SUM**2
                H_[8,5] = H_[5,8]
                H_[6,8] = E1*E7*E9/SUM**2
                H_[8,6] = H_[6,8]
                H_[7,8] = E1*E8*E9/SUM**2
                H_[8,7] = H_[7,8]
                H_[8,8] = (-1.0e+0+E9/SUM)*E1*E9/SUM
                H_[0,9] = (-1.0e+0+E1/SUM)*E1*E10/SUM
                H_[9,0] = H_[0,9]
                H_[1,9] = E1*E2*E10/SUM**2
                H_[9,1] = H_[1,9]
                H_[2,9] = E1*E3*E10/SUM**2
                H_[9,2] = H_[2,9]
                H_[3,9] = E1*E4*E10/SUM**2
                H_[9,3] = H_[3,9]
                H_[4,9] = E1*E5*E10/SUM**2
                H_[9,4] = H_[4,9]
                H_[5,9] = E1*E6*E10/SUM**2
                H_[9,5] = H_[5,9]
                H_[6,9] = E1*E7*E10/SUM**2
                H_[9,6] = H_[6,9]
                H_[7,9] = E1*E8*E10/SUM**2
                H_[9,7] = H_[7,9]
                H_[8,9] = E1*E9*E10/SUM**2
                H_[9,8] = H_[8,9]
                H_[9,9] = (-1.0e+0+E10/SUM)*E1*E10/SUM
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

