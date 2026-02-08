from s2mpjlib import *
class  HS116(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS116
#    *******
# 
#    A 3 stage membrane separation model
# 
#    Source: problem 116 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: B. Debarsy, Feb 1990.
# 
#    classification = "C-CLQR2-MN-13-14"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS116'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 13
        v_['1'] = 1
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
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7']])
        valA = np.append(valA,float(-0.002))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(0.002))
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10']])
        valA = np.append(valA,float(-1.262626))
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C6')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-0.03475))
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C7')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(-0.03475))
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C8')
        [ig,ig_,_] = s2mpj_ii('C9',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C9')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('C10',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C10')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(-500.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(500.0))
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C11')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C12')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(-0.03475))
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C13')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(-1.262626))
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C14')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(-1.262626))
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
        self.gconst = arrset(self.gconst,ig_['C3'],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['C4'],float(50.0))
        self.gconst = arrset(self.gconst,ig_['C9'],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['C11'],float(0.9))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        grange = arrset(grange,ig_['C4'],float(200.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.1)
        self.xupper = np.full((self.n,1),+float('inf'))
        self.xupper[ix_['X1']] = 1.0
        self.xupper[ix_['X2']] = 1.0
        self.xupper[ix_['X3']] = 1.0
        self.xlower[ix_['X4']] = 0.0001
        self.xupper[ix_['X4']] = 0.1
        self.xupper[ix_['X5']] = 0.9
        self.xupper[ix_['X6']] = 0.9
        self.xupper[ix_['X7']] = 1000.0
        self.xupper[ix_['X8']] = 1000.0
        self.xlower[ix_['X9']] = 500.0
        self.xupper[ix_['X9']] = 1000.0
        self.xupper[ix_['X10']] = 500.0
        self.xlower[ix_['X11']] = 1.0
        self.xupper[ix_['X11']] = 150.0
        self.xlower[ix_['X12']] = 0.0001
        self.xupper[ix_['X12']] = 150.0
        self.xlower[ix_['X13']] = 0.0001
        self.xupper[ix_['X13']] = 150.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(0.5)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(0.5)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(0.8)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(0.8))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(0.9)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(0.9)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(0.1)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(0.1))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(0.14)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(0.14)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(0.5)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X6'])[0],float(0.5))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(489.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X7']),float(489.0)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(80.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X8'])[0],float(80.0))
        if('X9' in ix_):
            self.x0[ix_['X9']] = float(650.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X9']),float(650.0)))
        if('X10' in ix_):
            self.x0[ix_['X10']] = float(450.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X10'])[0],float(450.0))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(150.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X11']),float(150.0)))
        if('X12' in ix_):
            self.x0[ix_['X12']] = float(150.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X12'])[0],float(150.0))
        if('X13' in ix_):
            self.x0[ix_['X13']] = float(150.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X13']),float(150.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.1),None,None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['C5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.231059))
        ig = ig_['C6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.975))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.00975))
        ig = ig_['C7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.975))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.00975))
        ig = ig_['C8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['C9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E13'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        ig = ig_['C10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E15'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E17'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['C11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E19'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        ig = ig_['C12']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E20'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.975))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E21'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.00975))
        ig = ig_['C13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E22'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.231059))
        ig = ig_['C14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.231059))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               97.588409
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nle+self.neq,self.m)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLQR2-MN-13-14"
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
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(self, nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

