from s2mpjlib import *
class  SIPOW3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SIPOW3
#    *********
# 
#    This is a discretization of a one sided approximation problem of
#    approximating the function xi * xi * eta by a linear polynomial
#    on the boundary of the unit square [0,1]x[0,1].
# 
#    Source: problem 3 in
#    M. J. D. Powell,
#    "Log barrier methods for semi-infinite programming calculations"
#    Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.
# 
#    SIF input: A. R. Conn and Nick Gould, August 1993
# 
#    classification = "C-CLLR2-AN-4-V"
# 
#    Problem variants: they are identified by the values of M (even)
# 
# IE M                   20 
# IE M                   100 
# IE M                   500 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SIPOW3'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 2000
        v_['1'] = 1
        v_['2'] = 2
        v_['RM'] = float(v_['M'])
        v_['RM/8'] = 0.125*v_['RM']
        v_['RM/4'] = 0.25*v_['RM']
        v_['3RM/8'] = 0.375*v_['RM']
        v_['M/8'] = int(np.fix(v_['RM/8']))
        v_['M/4'] = int(np.fix(v_['RM/4']))
        v_['3M/8'] = int(np.fix(v_['3RM/8']))
        v_['M/8+1'] = 1+v_['M/8']
        v_['M/4+1'] = 1+v_['M/4']
        v_['3M/8+1'] = 1+v_['3M/8']
        v_['M/2'] = int(np.fix(v_['M']/v_['2']))
        v_['M/2+1'] = 1+v_['M/2']
        v_['RM'] = float(v_['M'])
        v_['STEP'] = 8.0/v_['RM']
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['I'] = -1+J
            v_['RI'] = float(v_['I'])
            v_['XI'+str(J)] = v_['RI']*v_['STEP']
        for J in range(int(v_['1']),int(v_['M/8'])+1):
            v_['RJ'] = float(J)
            v_['ETA'+str(J)] = v_['XI'+str(J)]
            v_['XI'+str(J)] = 0.0
        for J in range(int(v_['M/8+1']),int(v_['M/4'])+1):
            v_['RJ'] = float(J)
            v_['XI'+str(J)] = -1.0+v_['XI'+str(J)]
            v_['ETA'+str(J)] = 1.0
        for J in range(int(v_['M/4+1']),int(v_['3M/8'])+1):
            v_['RJ'] = float(J)
            v_['ETA'+str(J)] = -2.0+v_['XI'+str(J)]
            v_['XI'+str(J)] = 1.0
        for J in range(int(v_['3M/8+1']),int(v_['M/2'])+1):
            v_['RJ'] = float(J)
            v_['XI'+str(J)] = -3.0+v_['XI'+str(J)]
            v_['ETA'+str(J)] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        self.xnames=arrset(self.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        self.xnames=arrset(self.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        self.xnames=arrset(self.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        self.xnames=arrset(self.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(1.0))
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(J))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X1']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X4']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X2']])
            valA = np.append(valA,float(v_['XI'+str(J)]))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X3']])
            valA = np.append(valA,float(v_['ETA'+str(J)]))
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['J+'] = v_['M/2']+J
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X1']])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X2']])
            valA = np.append(valA,float(v_['XI'+str(J)]))
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X3']])
            valA = np.append(valA,float(v_['ETA'+str(J)]))
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
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['J+'] = v_['M/2']+J
            v_['XIXI'] = v_['XI'+str(J)]*v_['XI'+str(J)]
            v_['XIXIETA'] = v_['XIXI']*v_['ETA'+str(J)]
            self.gconst = arrset(self.gconst,ig_['C'+str(J)],float(v_['XIXIETA']))
            self.gconst  = (
                  arrset(self.gconst,ig_['C'+str(int(v_['J+']))],float(v_['XIXIETA'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('X1' in ix_):
            self.x0[ix_['X1']] = float(-0.1)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X1']),float(-0.1)))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(0.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(0.0)))
        if('X4' in ix_):
            self.x0[ix_['X4']] = float(1.2)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X4']),float(1.2)))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            3.0315716D-1 ! m = 20
# LO SOLUTION            5.0397238D-1 ! m = 100
# LO SOLUTION            5.3016386D-1 ! m = 500
# LO SOLUTION            5.3465470D-1 ! m = 2000
# LO SOLUTION            5.3564207D-1 ! m = 10000
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CLLR2-AN-4-V"
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

