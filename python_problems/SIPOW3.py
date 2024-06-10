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
#    classification = "LLR2-AN-4-V"
# 
#    Problem variants: they are identified by the values of M (even)
# 
# IE M                   20 
# IE M                   100 
# IE M                   500 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SIPOW3'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(J))
            iv = ix_['X1']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X4']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X2']
            pbm.A[ig,iv] = float(v_['XI'+str(J)])+pbm.A[ig,iv]
            iv = ix_['X3']
            pbm.A[ig,iv] = float(v_['ETA'+str(J)])+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['J+'] = v_['M/2']+J
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            iv = ix_['X1']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            iv = ix_['X2']
            pbm.A[ig,iv] = float(v_['XI'+str(J)])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(int(v_['J+'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(int(v_['J+'])))
            iv = ix_['X3']
            pbm.A[ig,iv] = float(v_['ETA'+str(J)])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = find(gtype,lambda x:x=='<=')
        eqgrps = find(gtype,lambda x:x=='==')
        gegrps = find(gtype,lambda x:x=='>=')
        pb.nle = len(legrps)
        pb.neq = len(eqgrps)
        pb.nge = len(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = find(gtype,lambda x:(x=='<=' or x=='==' or x=='>='))
        pb.cnames= cnames[pbm.congrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = find(gtype,lambda x:x=='<>')
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['J+'] = v_['M/2']+J
            v_['XIXI'] = v_['XI'+str(J)]*v_['XI'+str(J)]
            v_['XIXIETA'] = v_['XIXI']*v_['ETA'+str(J)]
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(J)],float(v_['XIXIETA']))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['C'+str(int(v_['J+']))],float(v_['XIXIETA'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(-0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(-0.1)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(0.0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(0.0)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(1.2)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4']),float(1.2)))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            3.0315716D-1 ! m = 20
# LO SOLUTION            5.0397238D-1 ! m = 100
# LO SOLUTION            5.3016386D-1 ! m = 500
# LO SOLUTION            5.3465470D-1 ! m = 2000
# LO SOLUTION            5.3564207D-1 ! m = 10000
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "LLR2-AN-4-V"
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

