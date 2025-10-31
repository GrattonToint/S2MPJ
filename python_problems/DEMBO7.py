from s2mpjlib import *
class  DEMBO7(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEMBO7
#    *******
# 
#    A 7 stage membrane separation model
# 
#    Source: problem 7 in
#    R.S. Dembo,
#    "A set of geometric programming test problems and their solutions",
#    Mathematical Programming, 17, 192-213, 1976.
# 
#    SIF input: A. R. Conn, June 1993.
# 
#    classification = "C-CQOR2-MN-16-20"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DEMBO7'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 16
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
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(1.262626))
        [ig,ig_,_] = s2mpj_ii('C0',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C0')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X14']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X15']])
        valA = np.append(valA,float(1.262626))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(1.262626))
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(0.975))
        [ig,ig_,_] = s2mpj_ii('C2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(0.975))
        [ig,ig_,_] = s2mpj_ii('C3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(0.975))
        [ig,ig_,_] = s2mpj_ii('C4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C4')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(0.975))
        [ig,ig_,_] = s2mpj_ii('C5',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C5')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(0.975))
        [ig,ig_,_] = s2mpj_ii('C6',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C6')
        [ig,ig_,_] = s2mpj_ii('C7',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C7')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X13']])
        valA = np.append(valA,float(-0.002))
        [ig,ig_,_] = s2mpj_ii('C8',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C8')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('C9',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C9')
        [ig,ig_,_] = s2mpj_ii('C10',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C10')
        [ig,ig_,_] = s2mpj_ii('C11',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C11')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X16']])
        valA = np.append(valA,float(0.002))
        [ig,ig_,_] = s2mpj_ii('C12',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C12')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X11']])
        valA = np.append(valA,float(0.002))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X12']])
        valA = np.append(valA,float(-0.002))
        [ig,ig_,_] = s2mpj_ii('C13',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C13')
        [ig,ig_,_] = s2mpj_ii('C14',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C14')
        [ig,ig_,_] = s2mpj_ii('C15',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C15')
        [ig,ig_,_] = s2mpj_ii('C16',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C16')
        [ig,ig_,_] = s2mpj_ii('C17',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C17')
        [ig,ig_,_] = s2mpj_ii('C18',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C18')
        [ig,ig_,_] = s2mpj_ii('C19',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'C19')
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
        self.gconst = arrset(self.gconst,ig_['C0'],float(50.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = -np.full((self.nle,1),float('inf'))
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        grange = arrset(grange,ig_['C0'],float(200.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.1)
        self.xupper = np.full((self.n,1),0.9)
        self.xupper[ix_['X5']] = 1.0
        self.xlower[ix_['X5']] = 0.9
        self.xupper[ix_['X6']] = 0.1
        self.xlower[ix_['X6']] = 0.0001
        self.xupper[ix_['X11']] = 1000.0
        self.xlower[ix_['X11']] = 1.0
        self.xupper[ix_['X12']] = 500.0
        self.xlower[ix_['X12']] = 0.000001
        self.xupper[ix_['X13']] = 500.0
        self.xlower[ix_['X13']] = 1.0
        self.xupper[ix_['X14']] = 1000.0
        self.xlower[ix_['X14']] = 500.0
        self.xupper[ix_['X15']] = 1000.0
        self.xlower[ix_['X15']] = 500.0
        self.xupper[ix_['X16']] = 500.0
        self.xlower[ix_['X16']] = 0.000001
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(0.8)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(0.8)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(0.83)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X2'])[0],float(0.83))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(0.85)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(0.85)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(0.87)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X4'])[0],float(0.87))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(0.90)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(0.90)))
        if('X6' in ix_):
            self.x0[ix_['X6']] = float(0.10)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X6'])[0],float(0.10))
        if('X7' in ix_):
            self.x0[ix_['X7']] = float(0.12)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X7']),float(0.12)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(0.19)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X8'])[0],float(0.19))
        if('X9' in ix_):
            self.x0[ix_['X9']] = float(0.25)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X9']),float(0.25)))
        if('X10' in ix_):
            self.x0[ix_['X10']] = float(0.29)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X10'])[0],float(0.29))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(512.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X11']),float(512.0)))
        if('X12' in ix_):
            self.x0[ix_['X12']] = float(13.1)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X12'])[0],float(13.1))
        if('X13' in ix_):
            self.x0[ix_['X13']] = float(71.8)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X13']),float(71.8)))
        if('X14' in ix_):
            self.x0[ix_['X14']] = float(640.0)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X14'])[0],float(640.0))
        if('X15' in ix_):
            self.x0[ix_['X15']] = float(650.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X15']),float(650.0)))
        if('X16' in ix_):
            self.x0[ix_['X16']] = float(5.7)
        else:
            self.y0 = arrset(self.y0,np.where(self.congrps==ig_['X16'])[0],float(5.7))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eINV', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eQT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eQTQT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'W')
        elftv = loaset(elftv,it,3,'Z')
        [it,iet_,_] = s2mpj_ii( 'en2PRRC', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2mpj_ii( 'eQTRC', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2mpj_ii( 'eSQQT', iet_)
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
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQT')
        ielftype = arrset(ielftype,ie,iet_["eSQQT"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQT')
        ielftype = arrset(ielftype,ie,iet_["eSQQT"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQT')
        ielftype = arrset(ielftype,ie,iet_["eSQQT"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQT')
        ielftype = arrset(ielftype,ie,iet_["eSQQT"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E14'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E15'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQQT')
        ielftype = arrset(ielftype,ie,iet_["eSQQT"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E16'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTQT')
        ielftype = arrset(ielftype,ie,iet_["eQTQT"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTQT')
        ielftype = arrset(ielftype,ie,iet_["eQTQT"])
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PRRC')
        ielftype = arrset(ielftype,ie,iet_["en2PRRC"])
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E21'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PRRC')
        ielftype = arrset(ielftype,ie,iet_["en2PRRC"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E22'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PRRC')
        ielftype = arrset(ielftype,ie,iet_["en2PRRC"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E23'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E24'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E26'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PR')
        ielftype = arrset(ielftype,ie,iet_["en2PR"])
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E27'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E28'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTQT')
        ielftype = arrset(ielftype,ie,iet_["eQTQT"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E29'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTRC')
        ielftype = arrset(ielftype,ie,iet_["eQTRC"])
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E30'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTRC')
        ielftype = arrset(ielftype,ie,iet_["eQTRC"])
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E31'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTQT')
        ielftype = arrset(ielftype,ie,iet_["eQTQT"])
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E32'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTQT')
        ielftype = arrset(ielftype,ie,iet_["eQTQT"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E33'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E34'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eINV')
        ielftype = arrset(ielftype,ie,iet_["eINV"])
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E35'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E36'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQTRC')
        ielftype = arrset(ielftype,ie,iet_["eQTRC"])
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E37'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eINV')
        ielftype = arrset(ielftype,ie,iet_["eINV"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E38'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PRRC')
        ielftype = arrset(ielftype,ie,iet_["en2PRRC"])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E39'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E40'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E41'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E42'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E43'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E44'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E45'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQT')
        ielftype = arrset(ielftype,ie,iet_["eQT"])
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(0.1),float(0.9),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        ig = ig_['C0']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.231060))
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.034750))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.00975))
        ig = ig_['C2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.034750))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.00975))
        ig = ig_['C3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.034750))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.00975))
        ig = ig_['C4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E12'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.034750))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E13'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.00975))
        ig = ig_['C5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E14'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.034750))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E15'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.00975))
        ig = ig_['C6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E16'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E17'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E18'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E19'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E20'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E21'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E22'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        ig = ig_['C8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E24'])
        self.grelw = loaset(self.grelw,ig,posel,float(0.002))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E26'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        ig = ig_['C9']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E27'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E28'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E29'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(500.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E30'])
        self.grelw = loaset(self.grelw,ig,posel,float(-500.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E31'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C10']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E32'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E33'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E34'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(500.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E35'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E36'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-500.0))
        ig = ig_['C11']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E37'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.9))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E38'])
        self.grelw = loaset(self.grelw,ig,posel,float(-0.002))
        ig = ig_['C13']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E39'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C14']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E40'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C15']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E41'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C16']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E42'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C17']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E43'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C18']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E44'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        ig = ig_['C19']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E45'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               174.788807
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle)] = grange[legrps]
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nle+self.neq,self.m)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQOR2-MN-16-20"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eINV(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 1.0/EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/EV_[0]**2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0/EV_[0]**3
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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]/EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1/EV_[1]
            g_[1] = -EV_[0]/EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0/EV_[1]**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (2.0*EV_[0])/EV_[1]**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQTQT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XW = EV_[0]*EV_[2]
        YZ = EV_[1]*EV_[3]
        f_   = XW/YZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[2]/YZ
            g_[1] = -XW/(EV_[1]*YZ)
            g_[2] = EV_[0]/YZ
            g_[3] = -XW/(YZ*EV_[3])
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = -EV_[2]/(EV_[1]*YZ)
                H_[1,0] = H_[0,1]
                H_[0,2] = 1.0/YZ
                H_[2,0] = H_[0,2]
                H_[0,3] = -EV_[2]/(EV_[3]*YZ)
                H_[3,0] = H_[0,3]
                H_[1,1] = (2.0*XW)/(EV_[1]**2*YZ)
                H_[1,2] = -EV_[0]/(EV_[1]*YZ)
                H_[2,1] = H_[1,2]
                H_[1,3] = XW/YZ**2
                H_[3,1] = H_[1,3]
                H_[2,3] = -EV_[0]/(EV_[3]*YZ)
                H_[3,2] = H_[2,3]
                H_[3,3] = (2.0*XW)/(EV_[3]**2*YZ)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PRRC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]/EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]/EV_[2]
            g_[1] = EV_[0]/EV_[2]
            g_[2] = -EV_[0]*EV_[1]/EV_[2]**2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 1.0/EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = -EV_[1]/EV_[2]**2
                H_[2,0] = H_[0,2]
                H_[1,2] = -EV_[0]/EV_[2]**2
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[0]*EV_[1]/EV_[2]**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQTRC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        YZ = EV_[1]*EV_[2]
        f_   = EV_[0]/YZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/YZ
            g_[1] = -EV_[0]/(EV_[1]*YZ)
            g_[2] = -EV_[0]/(EV_[2]*YZ)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = -1.0/(EV_[1]*YZ)
                H_[1,0] = H_[0,1]
                H_[0,2] = -1.0/(EV_[2]*YZ)
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0*EV_[0])/(EV_[1]**2*YZ)
                H_[1,2] = EV_[0]/YZ**2
                H_[2,1] = H_[1,2]
                H_[2,2] = (2.0*EV_[0])/(EV_[2]**2*YZ)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQQT(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**2/EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]/EV_[1]
            g_[1] = -EV_[0]**2/EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0/EV_[1]
                H_[0,1] = -2.0*EV_[0]/EV_[1]**2
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0]**2/EV_[1]**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

