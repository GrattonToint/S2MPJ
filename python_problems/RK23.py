from s2mpjlib import *
class  RK23(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RK23
#    *********
# 
#    Find coefficients for an embedded pair of explicit 2nd
#    and 3rd order Runge Kutta Method such that the leading
#    term in the local truncation error is minimized.
# 
#    Source:
#    Similar ideas for 4th and 5th order pairs are discussed in:
#    Hairer, Norsett and Wanner, Solving Ordinary Differential
#    Equations I, Springer 1980, page 158 ff.
# 
#    SIF input: S. Leyffer, January 1997.
# 
#    classification = "C-CLOR2-RN-17-11"
# 
# 
#    ... COMPUTED PARAMETERS
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RK23'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['ONE'] = 1.0
        v_['THREE'] = 3.0
        v_['FOUR'] = 4.0
        v_['SIX'] = 6.0
        v_['ONETHIRD'] = v_['ONE']/v_['THREE']
        v_['ONESIXTH'] = v_['ONE']/v_['SIX']
        v_['FOURSIXTH'] = v_['FOUR']/v_['SIX']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('C2',ix_)
        self.xnames=arrset(self.xnames,iv,'C2')
        [iv,ix_,_] = s2mpj_ii('A21',ix_)
        self.xnames=arrset(self.xnames,iv,'A21')
        [iv,ix_,_] = s2mpj_ii('C3',ix_)
        self.xnames=arrset(self.xnames,iv,'C3')
        [iv,ix_,_] = s2mpj_ii('A31',ix_)
        self.xnames=arrset(self.xnames,iv,'A31')
        [iv,ix_,_] = s2mpj_ii('A32',ix_)
        self.xnames=arrset(self.xnames,iv,'A32')
        [iv,ix_,_] = s2mpj_ii('B1',ix_)
        self.xnames=arrset(self.xnames,iv,'B1')
        [iv,ix_,_] = s2mpj_ii('B2',ix_)
        self.xnames=arrset(self.xnames,iv,'B2')
        [iv,ix_,_] = s2mpj_ii('B3',ix_)
        self.xnames=arrset(self.xnames,iv,'B3')
        [iv,ix_,_] = s2mpj_ii('BB1',ix_)
        self.xnames=arrset(self.xnames,iv,'BB1')
        [iv,ix_,_] = s2mpj_ii('BB2',ix_)
        self.xnames=arrset(self.xnames,iv,'BB2')
        [iv,ix_,_] = s2mpj_ii('BB3',ix_)
        self.xnames=arrset(self.xnames,iv,'BB3')
        [iv,ix_,_] = s2mpj_ii('TP1',ix_)
        self.xnames=arrset(self.xnames,iv,'TP1')
        [iv,ix_,_] = s2mpj_ii('TM1',ix_)
        self.xnames=arrset(self.xnames,iv,'TM1')
        [iv,ix_,_] = s2mpj_ii('TP2',ix_)
        self.xnames=arrset(self.xnames,iv,'TP2')
        [iv,ix_,_] = s2mpj_ii('TM2',ix_)
        self.xnames=arrset(self.xnames,iv,'TM2')
        [iv,ix_,_] = s2mpj_ii('TP3',ix_)
        self.xnames=arrset(self.xnames,iv,'TP3')
        [iv,ix_,_] = s2mpj_ii('TM3',ix_)
        self.xnames=arrset(self.xnames,iv,'TM3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('ROWS1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ROWS1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['A21']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['C2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('ROWS2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ROWS2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['A31']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['A32']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['C3']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('FIRST2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'FIRST2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['B1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['B2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['B3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('FIRST3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'FIRST3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['BB1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['BB2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['BB3']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('SECND2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SECND2')
        [ig,ig_,_] = s2mpj_ii('SECND3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SECND3')
        [ig,ig_,_] = s2mpj_ii('THIRD31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'THIRD31')
        [ig,ig_,_] = s2mpj_ii('THIRD32',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'THIRD32')
        [ig,ig_,_] = s2mpj_ii('ART1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('ART2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP2']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('ART3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TP3']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['TM3']])
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['ROWS1'],float(0.0))
        self.gconst = arrset(self.gconst,ig_['ROWS2'],float(0.0))
        self.gconst = arrset(self.gconst,ig_['FIRST2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['FIRST3'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['SECND2'],float(0.5))
        self.gconst = arrset(self.gconst,ig_['SECND3'],float(0.5))
        self.gconst = arrset(self.gconst,ig_['THIRD31'],float(v_['ONETHIRD']))
        self.gconst = arrset(self.gconst,ig_['THIRD32'],float(v_['ONESIXTH']))
        self.gconst = arrset(self.gconst,ig_['ART1'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['ART2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['ART3'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['TP1']] = 0.0
        self.xlower[ix_['TM1']] = 0.0
        self.xlower[ix_['TP2']] = 0.0
        self.xlower[ix_['TM2']] = 0.0
        self.xlower[ix_['TP3']] = 0.0
        self.xlower[ix_['TM3']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('C2' in ix_):
            self.x0[ix_['C2']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['C2']),float(1.0)))
        if('A21' in ix_):
            self.x0[ix_['A21']] = float(1.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['A21']),float(1.0)))
        if('C3' in ix_):
            self.x0[ix_['C3']] = float(0.5)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['C3']),float(0.5)))
        if('A31' in ix_):
            self.x0[ix_['A31']] = float(0.25)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['A31']),float(0.25)))
        if('A32' in ix_):
            self.x0[ix_['A32']] = float(0.25)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['A32']),float(0.25)))
        if('B1' in ix_):
            self.x0[ix_['B1']] = float(0.5)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B1']),float(0.5)))
        if('B2' in ix_):
            self.x0[ix_['B2']] = float(0.5)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B2']),float(0.5)))
        if('B3' in ix_):
            self.x0[ix_['B3']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['B3']),float(0.0)))
        self.x0[ix_['BB1']] = float(v_['ONESIXTH'])
        self.x0[ix_['BB2']] = float(v_['ONESIXTH'])
        self.x0[ix_['BB3']] = float(v_['FOURSIXTH'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'ePRODS', iet_)
        elftv = loaset(elftv,it,0,'W1')
        elftv = loaset(elftv,it,1,'W2')
        [it,iet_,_] = s2mpj_ii( 'ePRODQ', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        [it,iet_,_] = s2mpj_ii( 'eTPROD', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        [it,iet_,_] = s2mpj_ii( 'eQPROD', iet_)
        elftv = loaset(elftv,it,0,'Z1')
        elftv = loaset(elftv,it,1,'Z2')
        elftv = loaset(elftv,it,2,'Z3')
        elftv = loaset(elftv,it,3,'Z4')
        [it,iet_,_] = s2mpj_ii( 'eTPRODS', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'B2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'B3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'BB2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD')
        ielftype = arrset(ielftype,ie,iet_["ePROD"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePRODS')
        ielftype = arrset(ielftype,ie,iet_["ePRODS"])
        vname = 'BB2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='W1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='W2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePRODS')
        ielftype = arrset(ielftype,ie,iet_["ePRODS"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='W1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='W2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eTPROD')
        ielftype = arrset(ielftype,ie,iet_["eTPROD"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Y1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Y2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Y3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePRODQ')
        ielftype = arrset(ielftype,ie,iet_["ePRODQ"])
        vname = 'BB2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePRODQ')
        ielftype = arrset(ielftype,ie,iet_["ePRODQ"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='X2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eQPROD')
        ielftype = arrset(ielftype,ie,iet_["eQPROD"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Z1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Z2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Z3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='Z4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eTPRODS')
        ielftype = arrset(ielftype,ie,iet_["eTPRODS"])
        vname = 'BB3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['SECND2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['SECND3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['THIRD31']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['THIRD32']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['ART1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(4.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E9'])
        self.grelw = loaset(self.grelw,ig,posel,float(4.0))
        ig = ig_['ART2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(8.0))
        ig = ig_['ART3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(12.0))
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
        self.pbclass   = "C-CLOR2-RN-17-11"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*(EV_[1]**2.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]**2
            g_[1] = EV_[0]*2.0*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 2.0*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*(EV_[1]**3.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]**3.0
            g_[1] = EV_[0]*3.0*(EV_[1]**2.0)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 3.0*(EV_[1]**2.0)
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0]*6.0*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTPROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQPROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]
            g_[1] = EV_[0]*EV_[2]*EV_[3]
            g_[2] = EV_[0]*EV_[1]*EV_[3]
            g_[3] = EV_[0]*EV_[1]*EV_[2]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2]*EV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0]*EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[0]*EV_[1]
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTPRODS(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*(EV_[2]**2.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*(EV_[2]**2.0)
            g_[1] = EV_[0]*(EV_[2]**2.0)
            g_[2] = EV_[0]*EV_[1]*2.0*EV_[2]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]**2.0
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*2.0*EV_[2]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*2.0*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[0]*EV_[1]*2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

