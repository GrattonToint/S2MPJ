from s2mpjlib import *
class  DALLASS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DALLASS
#    *********
# 
#    The small Dallas water distribution problem
#    The problem is also named "W30" in some references.
#    This is a nonlinear network problem with conditioning of
#    the order of 10**4.
# 
#    Source:
#    R. Dembo,
#    private communication, 1986.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "ONR2-MN-46-31"
# 
#    Number of arcs
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DALLASS'

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
        v_['N'] = 46
        v_['NODES'] = 31
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X42']
        pbm.A[ig,iv] = float(-6.38400e+02)+pbm.A[ig,iv]
        iv = ix_['X43']
        pbm.A[ig,iv] = float(-6.33000e+02)+pbm.A[ig,iv]
        iv = ix_['X44']
        pbm.A[ig,iv] = float(-5.54500e+02)+pbm.A[ig,iv]
        iv = ix_['X45']
        pbm.A[ig,iv] = float(-5.05000e+02)+pbm.A[ig,iv]
        iv = ix_['X46']
        pbm.A[ig,iv] = float(-4.36900e+02)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N1')
        iv = ix_['X46']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X41']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N2')
        iv = ix_['X45']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N3')
        iv = ix_['X44']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N4')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N5')
        iv = ix_['X16']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N6')
        iv = ix_['X7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N7')
        iv = ix_['X9']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N8')
        iv = ix_['X10']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N9')
        iv = ix_['X12']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N10')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N11')
        iv = ix_['X15']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N12')
        iv = ix_['X20']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N13')
        iv = ix_['X42']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N14')
        iv = ix_['X21']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N15')
        iv = ix_['X43']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X21']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N16')
        iv = ix_['X14']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X23']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X22']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N17')
        iv = ix_['X23']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X25']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X24']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N18')
        iv = ix_['X31']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X25']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X22']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X26']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N19')
        iv = ix_['X26']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X28']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X27']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N20')
        iv = ix_['X28']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N21')
        iv = ix_['X31']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X30']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X29']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N22')
        iv = ix_['X30']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X27']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N23')
        iv = ix_['X24']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X32']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N24',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N24')
        iv = ix_['X38']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X29']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X34']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X33']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N25')
        iv = ix_['X32']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X35']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N26',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N26')
        iv = ix_['X35']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X37']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X36']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N27')
        iv = ix_['X37']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X34']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N28',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N28')
        iv = ix_['X36']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X40']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X39']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X38']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N29',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N29')
        iv = ix_['X39']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X33']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N30',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N30')
        iv = ix_['X40']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X41']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('N31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'N31')
        iv = ix_['X46']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X45']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X44']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X43']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X42']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['N5'],float(2.80000))
        pbm.gconst = arrset(pbm.gconst,ig_['N7'],float(4.03000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N8'],float(5.92000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N9'],float(1.15600))
        pbm.gconst = arrset(pbm.gconst,ig_['N10'],float(2.00000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N11'],float(4.95000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N16'],float(3.13000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N17'],float(8.44000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N18'],float(3.31000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N19'],float(5.30000e-02))
        pbm.gconst = arrset(pbm.gconst,ig_['N21'],float(2.72000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N22'],float(8.83000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N23'],float(5.71000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N24'],float(7.55000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N26'],float(5.27000e-01))
        pbm.gconst = arrset(pbm.gconst,ig_['N29'],float(1.00000e-03))
        pbm.gconst = arrset(pbm.gconst,ig_['N31'],float(-1.01960e+01))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-2.00000e+02)
        pb.xupper = np.full((pb.n,1),2.00000e+02)
        pb.xlower[ix_['X1']] = 0.00000
        pb.xupper[ix_['X1']] = 2.11673e+01
        pb.xlower[ix_['X2']] = 0.00000
        pb.xupper[ix_['X2']] = 4.37635e+01
        pb.xlower[ix_['X3']] = 0.00000
        pb.xupper[ix_['X3']] = 3.28255e+01
        pb.xlower[ix_['X19']] = 0.00000
        pb.xupper[ix_['X19']] = 2.20120e+01
        pb.xlower[ix_['X21']] = 0.00000
        pb.xupper[ix_['X21']] = 1.36703e+01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(-2.00000e+02))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(2.11673e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(2.11673e+01)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(4.37635e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),float(4.37635e+01)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(3.28255e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(3.28255e+01)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(1.42109e-14)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),float(1.42109e-14)))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(1.68826e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(1.68826e+02)))
        if('X7' in ix_):
            pb.x0[ix_['X7']] = float(2.81745e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X7']),float(2.81745e+01)))
        if('X8' in ix_):
            pb.x0[ix_['X8']] = float(8.75603e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X8']),float(8.75603e+01)))
        if('X9' in ix_):
            pb.x0[ix_['X9']] = float(-5.93858e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X9']),float(-5.93858e+01)))
        if('X10' in ix_):
            pb.x0[ix_['X10']] = float(-5.97888e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X10']),float(-5.97888e+01)))
        if('X11' in ix_):
            pb.x0[ix_['X11']] = float(1.83383e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X11']),float(1.83383e+02)))
        if('X13' in ix_):
            pb.x0[ix_['X13']] = float(-1.68331e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X13']),float(-1.68331e+02)))
        if('X15' in ix_):
            pb.x0[ix_['X15']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X15']),float(2.00000e+02)))
        if('X16' in ix_):
            pb.x0[ix_['X16']] = float(2.00000e-01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X16']),float(2.00000e-01)))
        if('X17' in ix_):
            pb.x0[ix_['X17']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X17']),float(2.00000e+02)))
        if('X18' in ix_):
            pb.x0[ix_['X18']] = float(-7.67574e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X18']),float(-7.67574e+01)))
        if('X19' in ix_):
            pb.x0[ix_['X19']] = float(2.20120e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X19']),float(2.20120e+01)))
        if('X20' in ix_):
            pb.x0[ix_['X20']] = float(1.36703e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X20']),float(1.36703e+01)))
        if('X21' in ix_):
            pb.x0[ix_['X21']] = float(1.36703e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X21']),float(1.36703e+01)))
        if('X22' in ix_):
            pb.x0[ix_['X22']] = float(-1.98461e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X22']),float(-1.98461e+02)))
        if('X23' in ix_):
            pb.x0[ix_['X23']] = float(1.81531e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X23']),float(1.81531e+02)))
        if('X24' in ix_):
            pb.x0[ix_['X24']] = float(-1.93133e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X24']),float(-1.93133e+01)))
        if('X25' in ix_):
            pb.x0[ix_['X25']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X25']),float(2.00000e+02)))
        if('X26' in ix_):
            pb.x0[ix_['X26']] = float(-1.98792e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X26']),float(-1.98792e+02)))
        if('X27' in ix_):
            pb.x0[ix_['X27']] = float(1.15500)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X27']),float(1.15500)))
        if('X28' in ix_):
            pb.x0[ix_['X28']] = float(0.00000)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X28']),float(0.00000)))
        if('X29' in ix_):
            pb.x0[ix_['X29']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X29']),float(2.00000e+02)))
        if('X30' in ix_):
            pb.x0[ix_['X30']] = float(2.72000e-01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X30']),float(2.72000e-01)))
        if('X32' in ix_):
            pb.x0[ix_['X32']] = float(-1.98843e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X32']),float(-1.98843e+01)))
        if('X33' in ix_):
            pb.x0[ix_['X33']] = float(1.78834e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X33']),float(1.78834e+02)))
        if('X34' in ix_):
            pb.x0[ix_['X34']] = float(-1.79589e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X34']),float(-1.79589e+02)))
        if('X35' in ix_):
            pb.x0[ix_['X35']] = float(-1.98843e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X35']),float(-1.98843e+01)))
        if('X37' in ix_):
            pb.x0[ix_['X37']] = float(1.79589e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X37']),float(1.79589e+02)))
        if('X40' in ix_):
            pb.x0[ix_['X40']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X40']),float(2.00000e+02)))
        if('X41' in ix_):
            pb.x0[ix_['X41']] = float(2.00000e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X41']),float(2.00000e+02)))
        if('X42' in ix_):
            pb.x0[ix_['X42']] = float(9.87694e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X42']),float(9.87694e+01)))
        if('X43' in ix_):
            pb.x0[ix_['X43']] = float(1.36703e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X43']),float(1.36703e+01)))
        if('X44' in ix_):
            pb.x0[ix_['X44']] = float(3.28255e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X44']),float(3.28255e+01)))
        if('X45' in ix_):
            pb.x0[ix_['X45']] = float(4.37635e+01)
        else:
            pb.y0  = (
                  arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X45']),float(4.37635e+01)))
        if('X46' in ix_):
            pb.x0[ix_['X46']] = float(-1.78833e+02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X46']),float(-1.78833e+02)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eT1', iet_)
        elftv = loaset(elftv,it,0,'ARC')
        elftp = []
        elftp = loaset(elftp,it,0,'C1')
        elftp = loaset(elftp,it,1,'C2')
        elftp = loaset(elftp,it,2,'C3')
        [it,iet_,_] = s2mpj_ii( 'eT4', iet_)
        elftv = loaset(elftv,it,0,'ARC')
        elftp = loaset(elftp,it,0,'C1')
        elftp = loaset(elftp,it,1,'C2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'X1'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.48060e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.51200e+02))
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'X2'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.91526e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.46300e+01))
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'X3'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.07752e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.81400e+01))
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X4'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.90000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X5'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X6'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.63000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E7'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X7'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.10000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X8'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.45000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X9'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(7.40000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X10'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.00000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X11'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(8.00000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.07000e+02))
        ename = 'E12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X12'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.20000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.80000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X13'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.00000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.80000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X14'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.00000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X15'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.12200e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.30000e+02))
        ename = 'E16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X16'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.50000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.22000e+02))
        ename = 'E17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X17'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.10000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X18'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.80000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.18000e+02))
        ename = 'E19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'X19'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.84530e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.12970e+02))
        ename = 'E20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X20'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.60000e+04))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.80000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eT4')
        ielftype = arrset(ielftype, ie, iet_["eT4"])
        vname = 'X21'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.86880e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.60610e+02))
        ename = 'E22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X22'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.20000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.36100e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.30000e+02))
        ename = 'E23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X23'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.60000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X24'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X25'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.60000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X26'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.30000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X27'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.20000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.24000e+02))
        ename = 'E28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X28'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X29'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.90000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X30'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.80000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X31'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.70000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X32'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.10000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X33'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X34'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(4.30000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.13000e+02))
        ename = 'E35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X35'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.20000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(9.50000e+01))
        ename = 'E36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X36'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.80000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X37'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(2.40000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+02))
        ename = 'E38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X38'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.31000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        ename = 'E39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X39'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(6.65000e+02))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X40'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.10000e+03))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.60000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.20000e+02))
        ename = 'E41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eT1')
            ielftype = arrset( ielftype,ie,iet_['eT1'])
        vname = 'X41'
        [iv,ix_,pb]  = (
              s2mpj_nlx(vname,ix_,pb,1,-2.00000e+02,2.00000e+02,-2.00000e+02))
        posev = find(elftv[ielftype[ie]],lambda x:x=='ARC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='C1')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(3.23000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C2')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+01))
        posep = find(elftp[ielftype[ie]],lambda x:x=='C3')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.00000e+02))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E12'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E14'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E16'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E17'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E18'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E19'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E20'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E21'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E22'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E23'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E24'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E26'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E27'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E28'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E29'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E30'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E31'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E32'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E33'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E34'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E35'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E36'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E37'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E38'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E39'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E40'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E41'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.2393D+04
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "ONR2-MN-46-31"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eT1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = 850559.0e0/2.85*pbm.elpar[iel_][0]
        TMP = TMP/(pbm.elpar[iel_][2]**1.85)
        TMP = TMP/(pbm.elpar[iel_][1]**4.87)
        X = np.absolute(EV_[0])
        XEXP = X**0.85
        f_   = TMP*X**2*XEXP
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.85*TMP*EV_[0]*XEXP
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 5.2725*TMP*XEXP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eT4(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EPS2 = 1.0e-14
        SQC1 = np.sqrt(pbm.elpar[iel_][0])
        X = min(EV_[0],SQC1)
        TMP = pbm.elpar[iel_][1]*(pbm.elpar[iel_][0]-X*X)
        TMP = np.sqrt(max(TMP,EPS2))
        TMP2 = np.sqrt(pbm.elpar[iel_][1])*np.arcsin(X/SQC1)
        f_   = 0.5*(-X*TMP-pbm.elpar[iel_][0]*TMP2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -TMP
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = pbm.elpar[iel_][1]*X/TMP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

