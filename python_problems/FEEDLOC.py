from s2mpjlib import *
class  FEEDLOC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FEEDLOC
#    *********
# 
#    Feed tray location & determination of optimum number of trays 
#    in a distillation column
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-CLOR2-AN-90-259"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FEEDLOC'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['M'] = 2
        v_['NMAX'] = 12
        v_['NMAX-1'] = -1+v_['NMAX']
        v_['F'] = 100.0
        v_['AL1'] = 1.0
        v_['AL2'] = 5.13435
        v_['XF1'] = 0.80
        v_['XF2'] = 0.20
        v_['SPEC'] = 0.001
        v_['BIGM'] = 1000.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [iv,ix_,_] = s2mpj_ii('S'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'S'+str(I))
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'W'+str(I))
            [iv,ix_,_] = s2mpj_ii('Z'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Z'+str(I))
        [iv,ix_,_] = s2mpj_ii('N',ix_)
        self.xnames=arrset(self.xnames,iv,'N')
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'Y'+str(I)+','+str(J))
        [iv,ix_,_] = s2mpj_ii('L',ix_)
        self.xnames=arrset(self.xnames,iv,'L')
        [iv,ix_,_] = s2mpj_ii('V',ix_)
        self.xnames=arrset(self.xnames,iv,'V')
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        self.xnames=arrset(self.xnames,iv,'R')
        [iv,ix_,_] = s2mpj_ii('P1',ix_)
        self.xnames=arrset(self.xnames,iv,'P1')
        [iv,ix_,_] = s2mpj_ii('P2',ix_)
        self.xnames=arrset(self.xnames,iv,'P2')
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
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('FENTR',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'FENTR')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('NTRAY',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NTRAY')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('NDEF1',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NDEF1')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(1.0))
            v_['RI'] = float(I)
            [ig,ig_,_] = s2mpj_ii('NDEF2',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NDEF2')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(v_['RI']))
        [ig,ig_,_] = s2mpj_ii('NDEF1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NDEF1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['N']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('NDEF2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NDEF2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['N']])
        valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NMAX-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('NIL'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'NIL'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            v_['RI'] = float(I)
            [ig,ig_,_] = s2mpj_ii('ENTX',ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'ENTX')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['RI']))
            v_['RI'] = -1.0*v_['RI']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(v_['RI']))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('LASTX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'LASTX'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('ZNOT'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'ZNOT'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(1.0))
            for K in range(int(I),int(v_['NMAX'])+1):
                [ig,ig_,_] = s2mpj_ii('ZNOT'+str(I),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'ZNOT'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['S'+str(K)]])
                valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('FEEDX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'FEEDX'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('WNES1u'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'WNES1u'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('WNES2u'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'WNES2u'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('PE1'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE1'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('PE2'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE2'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('PE3'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE3'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('PE4'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE4'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('PE1'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE1'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('PE2'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE2'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('PE3'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE3'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('PE4'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE4'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(-1.0))
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('XNOT'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'XNOT'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Z'+str(I)]])
                valA = np.append(valA,float(-1.0))
                [ig,ig_,_] = s2mpj_ii('YNOT'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'YNOT'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Y'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['Z'+str(I)]])
                valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            v_['TEMP'] = -1.0*v_['AL1']
            [ig,ig_,_] = s2mpj_ii('PHEE'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PHEE'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['1']))]])
            valA = np.append(valA,float(v_['TEMP']))
        [ig,ig_,_] = s2mpj_ii('DEFL',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'DEFL')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['L']])
        valA = np.append(valA,float(1.0))
        for J in range(int(v_['1']),int(v_['M'])+1):
            v_['TEMP'] = -1.0*v_['F']
            [ig,ig_,_] = s2mpj_ii('CMB1u'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CMB1u'+str(J))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['2']))+','+str(J)]])
            valA = np.append(valA,float(v_['TEMP']))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                v_['TEMP'] = -1.0*v_['BIGM']
                [ig,ig_,_] = s2mpj_ii('CMBN1'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'CMBN1'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['S'+str(I)]])
                valA = np.append(valA,float(v_['BIGM']))
                [ig,ig_,_] = s2mpj_ii('CMBN2'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CMBN2'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['S'+str(I)]])
                valA = np.append(valA,float(v_['TEMP']))
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            v_['TEMP1'] = v_['F']*v_['XF'+str(int(v_['M']))]
            v_['TEMP1'] = -1.0*v_['TEMP1']
            [ig,ig_,_] = s2mpj_ii('CMB1'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CMB1'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(v_['TEMP']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(v_['BIGM']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['TEMP1']))
            [ig,ig_,_] = s2mpj_ii('CMB2'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CMB2'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(v_['BIGM']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['Z'+str(I)]])
            valA = np.append(valA,float(v_['TEMP']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(I)]])
            valA = np.append(valA,float(v_['TEMP1']))
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('RECR'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'RECR'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['S'+str(I)]])
            valA = np.append(valA,float(v_['BIGM']))
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
        self.gconst = arrset(self.gconst,ig_['FENTR'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['NTRAY'],float(1.0))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            self.gconst = arrset(self.gconst,ig_['WNES1u'+str(I)],float(1.0))
            self.gconst = arrset(self.gconst,ig_['WNES2u'+str(I)],float(1.0))
        v_['TEMP'] = -1.0*v_['BIGM']
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['CMBN1'+str(I)+','+str(J)],float(v_['BIGM'])))
                self.gconst  = (
                      arrset(self.gconst,ig_['CMBN2'+str(I)+','+str(J)],float(v_['TEMP'])))
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            self.gconst = arrset(self.gconst,ig_['CMB1'+str(I)],float(v_['BIGM']))
            v_['TEMP'] = -1.0*v_['BIGM']
            self.gconst = arrset(self.gconst,ig_['CMB2'+str(I)],float(v_['TEMP']))
        v_['TEMP'] = v_['XF'+str(int(v_['1']))]*v_['SPEC']
        v_['TEMP1'] = v_['TEMP']*v_['F']
        v_['RHS'] = v_['TEMP1']+v_['BIGM']
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            self.gconst = arrset(self.gconst,ig_['RECR'+str(I)],float(v_['RHS']))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = -np.full((self.nle,1),float('inf'))
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            grange = arrset(grange,ig_['PE1'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE2'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE3'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE4'+str(I)],float(2.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            self.xupper[ix_['Z'+str(I)]] = 1.0
            self.xupper[ix_['W'+str(I)]] = 1.0
            self.xupper[ix_['S'+str(I)]] = 1.0
            for J in range(int(v_['1']),int(v_['M'])+1):
                self.xupper[ix_['X'+str(I)+','+str(J)]] = 1.0
                self.xupper[ix_['Y'+str(I)+','+str(J)]] = 1.0
        self.xlower[ix_['N']] = 3.0
        v_['TEMP'] = float(v_['NMAX'])
        self.xupper[ix_['N']] = v_['TEMP']
        self.xlower[ix_['P2']] = 80.0
        self.xupper[ix_['P2']] = 80.0
        self.xupper[ix_['L']] = v_['F']
        self.xupper[ix_['V']] = v_['F']
        self.xupper[ix_['P1']] = v_['F']
        self.xupper[ix_['R']] = 5.0
        self.xlower[ix_['W'+str(int(v_['1']))]] = 0.0
        self.xupper[ix_['W'+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['W'+str(int(v_['2']))]] = 0.0
        self.xupper[ix_['W'+str(int(v_['2']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.5))
        self.y0 = np.full((self.m,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eA2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                ename = 'PHE'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['AL'+str(K)]))
        ename = 'DEFLE'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eA2PROD')
        ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
        vname = 'R'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'P1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='A')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'CMB11u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'P2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            ename = 'CMB12u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['1']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            ename = 'CMB13u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['2']))+','+str(J)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['M'])+1):
                ename = 'CM11'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'L'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                ename = 'CM12'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'P1'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                ename = 'CM13'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'V'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
                ename = 'CM21'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'L'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                ename = 'CM22'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'P1'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
                ename = 'CM23'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'V'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ename = 'C11'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            for K in range(int(I),int(v_['NMAX'])+1):
                ename = 'C12'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['F']))
            ename = 'C13'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            ename = 'C14'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                v_['TEMP'] = -1.0*v_['F']
                ename = 'C15'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TEMP']))
            ename = 'C16'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
            ename = 'C21'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            for K in range(int(v_['1']),int(v_['NMAX'])+1):
                ename = 'C22'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['F']))
            ename = 'C23'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
            ename = 'C24'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                v_['TEMP'] = -1.0*v_['F']
                ename = 'C25'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='A')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['TEMP']))
            ename = 'C26'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            ename = 'REC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype,ie,iet_["eA2PROD"])
            vname = 'P1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.5))
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='A')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['PHEE'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['PHE'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['DEFL']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['DEFLE'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['CMB1u'+str(J)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CMB11u'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['CMB12u'+str(J)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['CMB13u'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['CMBN1'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM11'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM12'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM13'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['CMBN2'+str(I)+','+str(J)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM21'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                posel = posel+1
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM22'+str(I)+','+str(J)])
                self.grelw = loaset(self.grelw,ig,posel, 1.)
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['CM23'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            ig = ig_['CMB1'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C11'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['C13'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C14'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['C16'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['CMB2'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C21'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['C23'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['C24'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['C26'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            for K in range(int(I),int(v_['NMAX'])+1):
                ig = ig_['CMB1'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C12'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['CMB2'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C22'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
            v_['I+1'] = 1+I
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                ig = ig_['CMB1'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C15'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
                ig = ig_['CMB2'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['C25'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            ig = ig_['RECR'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['REC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle)] = grange[legrps]
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nle+self.neq,self.m)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-AN-90-259"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eA2PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*EV_[1]
            g_[1] = self.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = self.elpar[iel_][0]
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

