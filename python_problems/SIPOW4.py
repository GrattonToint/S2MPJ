from s2mpjlib import *
class  SIPOW4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SIPOW4
#    *********
# 
#    This is a discretization of a one sided approximation problem of
#    approximating the function xi * xi * eta by a linear polynomial
#    on the boundary of a circle (xi - 0.5)**2 + (eta - 0.5)**2 = 0.5
# 
#    Source: problem 4 in
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

    name = 'SIPOW4'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'SIPOW4'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 2000
        v_['1'] = 1
        v_['2'] = 2
        v_['M/2'] = int(np.fix(v_['M']/v_['2']))
        v_['M/2+1'] = 1+v_['M/2']
        v_['RM'] = float(v_['M'])
        v_['1/RM'] = 1.0/v_['RM']
        v_['ONE'] = 1.0
        v_['HALF'] = 0.5
        v_['ROOTHALF'] = np.sqrt(v_['HALF'])
        v_['PI/4'] = np.arctan(v_['ONE'])
        v_['2PI'] = 8.0*v_['PI/4']
        v_['2PI/M'] = v_['2PI']*v_['1/RM']
        for J in range(int(v_['1']),int(v_['M/2'])+1):
            v_['RJ'] = float(J)
            v_['THETA'] = v_['RJ']*v_['2PI/M']
            v_['PI/4-T'] = v_['PI/4']-v_['THETA']
            v_['COS'] = np.cos(v_['PI/4-T'])
            v_['SIN'] = np.sin(v_['PI/4-T'])
            v_['RTC'] = v_['COS']*v_['ROOTHALF']
            v_['RTS'] = v_['SIN']*v_['ROOTHALF']
            v_['-RTC'] = -1.0*v_['RTC']
            v_['-RTS'] = -1.0*v_['RTS']
            v_['XI'+str(J)] = v_['HALF']+v_['-RTC']
            v_['ETA'+str(J)] = v_['HALF']+v_['-RTS']
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
# LO SOLUTION            2.0704432D-1 ! m = 20
# LO SOLUTION            2.6110334D-1 ! m = 100
# LO SOLUTION            2.7060094D-1 ! m = 500
# LO SOLUTION            2.7236200D-1 ! m = 2000
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

