from s2mpjlib import *
class  LOADBAL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The problem arises in the field of computer networks and parallel
#    computation.  It deals with the static load balancing in a tree
#    computer network with two-way traffic.  A set of heterogeneous host
#    computers are interconnected, in which each node processes jobs (the 
#    jobs arriving at each node according to a time invariant Poisson process) 
#    locally or sends it to a remote node,.  In the latter case, there is a
#    communication delay of forwarding the job and getting a response back.
#    The problem is then to minimize the mean response time of a job.
# 
#    The example considered here features 11 computers arranged as follows:
# 
#          1      6      9
#           \     |     /
#            \    |    /
#         2---4---5---8---10
#            /    |    \
#           /     |     \
#          3      7      11
# 
#    Source:
#    J. Li and H. Kameda,
#    "Optimal load balancing in tree network with two-way traffic",
#    Computer networks and ISDN systems, vol. 25, pp. 1335-1348, 1993.
# 
#    SIF input: Masha Sosonkina, Virginia Tech., 1995.
# 
#    classification = "OLR2-MN-31-31"
# 
#  Parameter assignment.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LOADBAL'

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
        v_['1'] = 1
        v_['P1'] = 3
        v_['N'] = 11
        v_['NLINK'] = 20
        v_['NLINK-3'] = 17
        v_['NLINK-4'] = 16
        v_['4C'] = 4
        v_['5C'] = 5
        v_['6C'] = 6
        v_['7C'] = 7
        v_['8C'] = 8
        v_['FI'] = 514.0
        v_['0.2*FI'] = 0.2*v_['FI']
        v_['0.0125*FI'] = 0.0125*v_['FI']
        v_['0.05*FI'] = 0.05*v_['FI']
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['0.2*FI']))
        for I in range(int(v_['1']),int(v_['NLINK'])+1):
            [ig,ig_,_] = s2mpj_ii('CNST'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CNST'+str(I))
            [ig,ig_,_] = s2mpj_ii('GA'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['0.0125*FI']))
            [ig,ig_,_] = s2mpj_ii('GB'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['0.05*FI']))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('N'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'N'+str(I))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('X4,1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,1')
        ig = ig_['N1']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X1,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1,4')
        ig = ig_['N1']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,2')
        ig = ig_['N2']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X2,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2,4')
        ig = ig_['N2']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,3')
        ig = ig_['N3']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X3,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3,4')
        ig = ig_['N3']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,5')
        ig = ig_['N5']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,4')
        ig = ig_['N5']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N4']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,6',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,6')
        ig = ig_['N6']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X6,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X6,5')
        ig = ig_['N6']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,7',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,7')
        ig = ig_['N7']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X7,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X7,5')
        ig = ig_['N7']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,8')
        ig = ig_['N8']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,5')
        ig = ig_['N8']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N5']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,9',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,9')
        ig = ig_['N9']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X9,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X9,8')
        ig = ig_['N9']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,10',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,10')
        ig = ig_['N10']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X10,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X10,8')
        ig = ig_['N10']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,11',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,11')
        ig = ig_['N11']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X11,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X11,8')
        ig = ig_['N11']
        pbm.A[ig,iv] = float(-1.00)+pbm.A[ig,iv]
        ig = ig_['N8']
        pbm.A[ig,iv] = float(1.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,1')
        ig = ig_['CNST1']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST2']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X1,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1,4')
        ig = ig_['CNST1']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST2']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,2')
        ig = ig_['CNST3']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST4']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X2,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2,4')
        ig = ig_['CNST3']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST4']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,3')
        ig = ig_['CNST5']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST6']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X3,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3,4')
        ig = ig_['CNST5']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST6']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,6',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,6')
        ig = ig_['CNST7']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST8']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X6,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X6,5')
        ig = ig_['CNST7']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST8']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,7',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,7')
        ig = ig_['CNST9']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST10']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X7,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X7,5')
        ig = ig_['CNST9']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST10']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,9',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,9')
        ig = ig_['CNST11']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST12']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X9,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X9,8')
        ig = ig_['CNST11']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST12']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,10',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,10')
        ig = ig_['CNST13']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST14']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X10,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X10,8')
        ig = ig_['CNST13']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        ig = ig_['CNST14']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,11',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,11')
        ig = ig_['CNST15']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        ig = ig_['CNST16']
        pbm.A[ig,iv] = float(80.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X11,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X11,8')
        ig = ig_['CNST15']
        pbm.A[ig,iv] = float(80.0)+pbm.A[ig,iv]
        ig = ig_['CNST16']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X4,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4,5')
        ig = ig_['CNST17']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        ig = ig_['CNST18']
        pbm.A[ig,iv] = float(80.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,4')
        ig = ig_['CNST17']
        pbm.A[ig,iv] = float(80.0)+pbm.A[ig,iv]
        ig = ig_['CNST18']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X5,8',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5,8')
        ig = ig_['CNST19']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        ig = ig_['CNST20']
        pbm.A[ig,iv] = float(80.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X8,5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X8,5')
        ig = ig_['CNST19']
        pbm.A[ig,iv] = float(80.0)+pbm.A[ig,iv]
        ig = ig_['CNST20']
        pbm.A[ig,iv] = float(20.00)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
            ig = ig_['N'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
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
        pbm.gconst = arrset(pbm.gconst,ig_['N1'],float(-95.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N2'],float(-95.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N3'],float(-19.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N4'],float(-70.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N5'],float(-70.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N6'],float(-19.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N7'],float(-19.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N8'],float(-70.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N9'],float(-19.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N10'],float(-19.0))
        pbm.gconst = arrset(pbm.gconst,ig_['N11'],float(-19.0))
        v_['CIJE'] = 999.99
        for I in range(int(v_['1']),int(v_['NLINK-4'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CNST'+str(I)],float(v_['CIJE']))
        v_['CIJE'] = 9999.99
        for I in range(int(v_['NLINK-3']),int(v_['NLINK'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CNST'+str(I)],float(v_['CIJE']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['B1']] = 99.99
        pb.xupper[ix_['B2']] = 99.99
        pb.xupper[ix_['B4']] = 99.99
        pb.xupper[ix_['B5']] = 99.99
        pb.xupper[ix_['B8']] = 99.99
        pb.xupper[ix_['B3']] = 19.99
        pb.xupper[ix_['B6']] = 19.99
        pb.xupper[ix_['B7']] = 19.99
        pb.xupper[ix_['B9']] = 19.99
        pb.xupper[ix_['B10']] = 19.99
        pb.xupper[ix_['B11']] = 19.99
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1,4' in ix_):
            pb.x0[ix_['X1,4']] = float(00.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1,4']),float(00.0)))
        if('X4,1' in ix_):
            pb.x0[ix_['X4,1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,1']),float(0.0)))
        if('B1' in ix_):
            pb.x0[ix_['B1']] = float(95.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B1']),float(95.0)))
        if('X2,4' in ix_):
            pb.x0[ix_['X2,4']] = float(0.00)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2,4']),float(0.00)))
        if('X4,2' in ix_):
            pb.x0[ix_['X4,2']] = float(0.00)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,2']),float(0.00)))
        if('B2' in ix_):
            pb.x0[ix_['B2']] = float(95.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B2']),float(95.0)))
        if('X3,4' in ix_):
            pb.x0[ix_['X3,4']] = float(0.00)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3,4']),float(0.00)))
        if('X4,3' in ix_):
            pb.x0[ix_['X4,3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,3']),float(0.0)))
        if('B3' in ix_):
            pb.x0[ix_['B3']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B3']),float(19.0)))
        if('B4' in ix_):
            pb.x0[ix_['B4']] = float(70.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B4']),float(70.0)))
        if('X5,4' in ix_):
            pb.x0[ix_['X5,4']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,4']),float(0.0)))
        if('X4,5' in ix_):
            pb.x0[ix_['X4,5']] = float(00.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,5']),float(00.0)))
        if('X6,5' in ix_):
            pb.x0[ix_['X6,5']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X6,5']),float(0.0)))
        if('X5,6' in ix_):
            pb.x0[ix_['X5,6']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,6']),float(0.0)))
        if('X7,5' in ix_):
            pb.x0[ix_['X7,5']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X7,5']),float(0.0)))
        if('X5,7' in ix_):
            pb.x0[ix_['X5,7']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,7']),float(0.0)))
        if('B5' in ix_):
            pb.x0[ix_['B5']] = float(70.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B5']),float(70.0)))
        if('B6' in ix_):
            pb.x0[ix_['B6']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B6']),float(19.0)))
        if('B7' in ix_):
            pb.x0[ix_['B7']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B7']),float(19.0)))
        if('X8,5' in ix_):
            pb.x0[ix_['X8,5']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8,5']),float(0.0)))
        if('X5,8' in ix_):
            pb.x0[ix_['X5,8']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,8']),float(0.0)))
        if('X9,8' in ix_):
            pb.x0[ix_['X9,8']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X9,8']),float(0.0)))
        if('X8,9' in ix_):
            pb.x0[ix_['X8,9']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8,9']),float(0.0)))
        if('X10,8' in ix_):
            pb.x0[ix_['X10,8']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X10,8']),float(0.0)))
        if('X8,10' in ix_):
            pb.x0[ix_['X8,10']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8,10']),float(0.0)))
        if('X11,8' in ix_):
            pb.x0[ix_['X11,8']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X11,8']),float(0.0)))
        if('X8,11' in ix_):
            pb.x0[ix_['X8,11']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8,11']),float(0.0)))
        if('B8' in ix_):
            pb.x0[ix_['B8']] = float(70.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B8']),float(70.0)))
        if('B9' in ix_):
            pb.x0[ix_['B9']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B9']),float(19.0)))
        if('B10' in ix_):
            pb.x0[ix_['B10']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B10']),float(19.0)))
        if('B11' in ix_):
            pb.x0[ix_['B11']] = float(19.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B11']),float(19.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eBETA1', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eBETA2', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eCOMA1', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        [it,iet_,_] = s2mpj_ii( 'eCOMA2', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        [it,iet_,_] = s2mpj_ii( 'eCOMB1', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        [it,iet_,_] = s2mpj_ii( 'eCOMB2', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'EB1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA1')
        ielftype = arrset(ielftype, ie, iet_["eBETA1"])
        vname = 'B1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA1')
        ielftype = arrset(ielftype, ie, iet_["eBETA1"])
        vname = 'B2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA1')
        ielftype = arrset(ielftype, ie, iet_["eBETA1"])
        vname = 'B4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA1')
        ielftype = arrset(ielftype, ie, iet_["eBETA1"])
        vname = 'B5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA1')
        ielftype = arrset(ielftype, ie, iet_["eBETA1"])
        vname = 'B8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EB11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBETA2')
        ielftype = arrset(ielftype, ie, iet_["eBETA2"])
        vname = 'B11'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['P1'])+1):
            ename = 'EGA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            vname = 'X'+str(int(v_['4C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['4C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            vname = 'X'+str(int(v_['4C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['4C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['I+3'] = 3+I
            ename = 'EGA'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            ename = 'EGA'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['4C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGA'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['4C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            ename = 'EGB'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['4C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['4C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['I+6'] = 6+I
            v_['I+8'] = 8+I
            ename = 'EGA'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            ename = 'EGA'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['8C']))+','+str(int(v_['I+8']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGA'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+8']))+','+str(int(v_['8C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            ename = 'EGB'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['8C']))+','+str(int(v_['I+8']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+6']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+8']))+','+str(int(v_['8C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['I+9'] = 9+I
            ename = 'EGA'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            ename = 'EGA'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+8']))+','+str(int(v_['8C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGA'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['8C']))+','+str(int(v_['I+8']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            ename = 'EGB'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['I+8']))+','+str(int(v_['8C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I+9']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['8C']))+','+str(int(v_['I+8']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['6C']),int(v_['7C'])+1):
            v_['I2'] = 2*I
            v_['I2+1'] = 1+v_['I2']
            ename = 'EGA'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            ename = 'EGA'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['5C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGA'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['5C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            ename = 'EGB'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['5C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I2+1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['5C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['I2+2'] = 2+v_['I2']
            ename = 'EGA'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMA1')
            ielftype = arrset(ielftype, ie, iet_["eCOMA1"])
            ename = 'EGA'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['5C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGA'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['5C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCOMB1')
            ielftype = arrset(ielftype, ie, iet_["eCOMB1"])
            ename = 'EGB'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(I)+','+str(int(v_['5C']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EGB'+str(int(v_['I2+2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X'+str(int(v_['5C']))+','+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGA17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMA2')
        ielftype = arrset(ielftype, ie, iet_["eCOMA2"])
        vname = 'X5,4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGB17'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMB2')
        ielftype = arrset(ielftype, ie, iet_["eCOMB2"])
        vname = 'X5,4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGA18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMA2')
        ielftype = arrset(ielftype, ie, iet_["eCOMA2"])
        vname = 'X4,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5,4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGB18'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMB2')
        ielftype = arrset(ielftype, ie, iet_["eCOMB2"])
        vname = 'X4,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5,4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGA19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMA2')
        ielftype = arrset(ielftype, ie, iet_["eCOMA2"])
        vname = 'X5,8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGB19'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMB2')
        ielftype = arrset(ielftype, ie, iet_["eCOMB2"])
        vname = 'X5,8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGA20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMA2')
        ielftype = arrset(ielftype, ie, iet_["eCOMA2"])
        vname = 'X8,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5,8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EGB20'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCOMB2')
        ielftype = arrset(ielftype, ie, iet_["eCOMB2"])
        vname = 'X8,5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5,8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NLINK'])+1):
            ig = ig_['GB'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EGB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['GA'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EGA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-MN-31-31"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE ELEMENTS *
#  ROUTINE             *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,80.0)
        pbm.efpar = arrset( pbm.efpar,1,20.0)
        return pbm

    @staticmethod
    def eBETA2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CB = 20.0
        f_   = EV_[0]/(CB-EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = CB/((CB-EV_[0])**2)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2*CB/((CB-EV_[0])**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eBETA1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CB = 100.0
        f_   = EV_[0]/(CB-EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = CB/((CB-EV_[0])**2)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2*CB/((CB-EV_[0])**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOMA1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CIJ = 1000.0
        f_   = EV_[0]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  (CIJ-pbm.efpar[1]*EV_[1])/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**2)
            g_[1]  = (
                  EV_[0]*pbm.efpar[1]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      2*pbm.efpar[0]*(CIJ-pbm.efpar[1]*EV_[1])/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
                H_[0,1] = (pbm.efpar[1]*(CIJ+pbm.efpar[0]*EV_[0]-pbm.efpar[1]*EV_[1])/
                     (CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      2*pbm.efpar[1]*pbm.efpar[1]*EV_[0]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOMB1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CIJ = 1000.0
        f_   = EV_[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  (CIJ-pbm.efpar[0]*EV_[1])/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**2)
            g_[1]  = (
                  EV_[0]*pbm.efpar[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      2*pbm.efpar[1]*(CIJ-pbm.efpar[0]*EV_[1])/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
                H_[0,1] = (pbm.efpar[0]*(CIJ+pbm.efpar[1]*EV_[0]-pbm.efpar[0]*EV_[1])/
                     (CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      2*pbm.efpar[0]*pbm.efpar[0]*EV_[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOMA2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CIJ = 10000.0
        f_   = EV_[0]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  (CIJ-pbm.efpar[1]*EV_[1])/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**2)
            g_[1]  = (
                  EV_[0]*pbm.efpar[1]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      2*pbm.efpar[0]*(CIJ-pbm.efpar[1]*EV_[1])/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
                H_[0,1] = (pbm.efpar[1]*(CIJ+pbm.efpar[0]*EV_[0]-pbm.efpar[1]*EV_[1])/
                     (CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      2*pbm.efpar[1]*pbm.efpar[1]*EV_[0]/(CIJ-(pbm.efpar[0]*EV_[0]+pbm.efpar[1]*EV_[1]))**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCOMB2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        CIJ = 10000.0
        f_   = EV_[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  (CIJ-pbm.efpar[0]*EV_[1])/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**2)
            g_[1]  = (
                  EV_[0]*pbm.efpar[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0]  = (
                      2*pbm.efpar[1]*(CIJ-pbm.efpar[0]*EV_[1])/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
                H_[0,1] = (pbm.efpar[0]*(CIJ+pbm.efpar[1]*EV_[0]-pbm.efpar[0]*EV_[1])/
                     (CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
                H_[1,0] = H_[0,1]
                H_[1,1]  = (
                      2*pbm.efpar[0]*pbm.efpar[0]*EV_[0]/(CIJ-(pbm.efpar[1]*EV_[0]+pbm.efpar[0]*EV_[1]))**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

