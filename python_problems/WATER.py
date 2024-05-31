from s2xlib import *
class  WATER(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A small nonlinear network problem.
#    The problem is to compute the flows in a water distribution network
#    with 7 nodes and 8 links, subject to known supply/demand at the nodes 
#    and a unique reservoir at node 1.
# 
#    The problem is convex.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: E. P. Smith, Virginia Tech., Spring 1993.
# 
#    classification = "ONR2-MN-31-10"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'WATER'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'WATER'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('obj0102',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(105665.6))
        [ig,ig_,_] = s2x_ii('obj0203',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(3613.412))
        [ig,ig_,_] = s2x_ii('obj0204',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(105665.6))
        [ig,ig_,_] = s2x_ii('obj0305',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(890.1553))
        [ig,ig_,_] = s2x_ii('obj0405',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(76.66088))
        [ig,ig_,_] = s2x_ii('obj0406',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(55145.82))
        [ig,ig_,_] = s2x_ii('obj0607',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(26030.46))
        [ig,ig_,_] = s2x_ii('obj0705',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(890.1553))
        [ig,ig_,_] = s2x_ii('obj',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('c1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c1')
        [ig,ig_,_] = s2x_ii('c2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c2')
        [ig,ig_,_] = s2x_ii('c3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c3')
        [ig,ig_,_] = s2x_ii('c4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c4')
        [ig,ig_,_] = s2x_ii('c5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c5')
        [ig,ig_,_] = s2x_ii('c6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c6')
        [ig,ig_,_] = s2x_ii('c7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c7')
        [ig,ig_,_] = s2x_ii('c8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c8')
        [ig,ig_,_] = s2x_ii('c9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c9')
        [ig,ig_,_] = s2x_ii('c10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c10')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2x_ii('Q0102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0102')
        ig = ig_['obj0102']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c1']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0102')
        ig = ig_['c2']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0203')
        ig = ig_['obj0203']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c2']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0203')
        ig = ig_['c3']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0204')
        ig = ig_['obj0204']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c2']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0204')
        ig = ig_['c4']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0305',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0305')
        ig = ig_['obj0305']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c3']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0305',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0305')
        ig = ig_['c5']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0405',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0405')
        ig = ig_['obj0405']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c4']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0405',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0405')
        ig = ig_['c5']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0406',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0406')
        ig = ig_['obj0406']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c4']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0406',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0406')
        ig = ig_['c6']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0607',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0607')
        ig = ig_['obj0607']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c6']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0607',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0607')
        ig = ig_['c7']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0705',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0705')
        ig = ig_['obj0705']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c5']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0705',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0705')
        ig = ig_['c7']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q01u0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q01u0')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c1']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q01u0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q01u0')
        ig = ig_['c8']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y02up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y02up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c2']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y02up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y02up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y03up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y03up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c3']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y03up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y03up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y04up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y04up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c4']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y04up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y04up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y05up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y05up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c5']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y05up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y05up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y06up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y06up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c6']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y06up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y06up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y07up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y07up')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(210)+pbm.A[ig,iv]
        ig = ig_['c7']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('y07up',ix_)
        pb.xnames=arrset(pb.xnames,iv,'y07up')
        ig = ig_['c9']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu02',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu02')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-175)+pbm.A[ig,iv]
        ig = ig_['c2']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu02',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu02')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu03',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu03')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-190)+pbm.A[ig,iv]
        ig = ig_['c3']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu03',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu03')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu04',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu04')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-185)+pbm.A[ig,iv]
        ig = ig_['c4']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu04',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu04')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu05',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu05')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-180)+pbm.A[ig,iv]
        ig = ig_['c5']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu05',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu05')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu06',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu06')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-195)+pbm.A[ig,iv]
        ig = ig_['c6']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu06',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu06')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu07',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu07')
        ig = ig_['obj']
        pbm.A[ig,iv] = float(-190)+pbm.A[ig,iv]
        ig = ig_['c7']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yqu07',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yqu07')
        ig = ig_['c10']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0201',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0201')
        ig = ig_['c1']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c2']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0302',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0302')
        ig = ig_['c2']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c3']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0402',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0402')
        ig = ig_['c2']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c4']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0503',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0503')
        ig = ig_['c3']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c5']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0504',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0504')
        ig = ig_['c4']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c5']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0604',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0604')
        ig = ig_['c4']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c6']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0507',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0507')
        ig = ig_['c5']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c7']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Q0706',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q0706')
        ig = ig_['c6']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c7']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yupu0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yupu0')
        ig = ig_['c8']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        ig = ig_['c9']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('yu0uq',ix_)
        pb.xnames=arrset(pb.xnames,iv,'yu0uq')
        ig = ig_['c8']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['c10']
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['c1'],float(1120))
        pbm.gconst = arrset(pbm.gconst,ig_['c2'],float(-100))
        pbm.gconst = arrset(pbm.gconst,ig_['c3'],float(-100))
        pbm.gconst = arrset(pbm.gconst,ig_['c4'],float(-120))
        pbm.gconst = arrset(pbm.gconst,ig_['c5'],float(-270))
        pbm.gconst = arrset(pbm.gconst,ig_['c6'],float(-330))
        pbm.gconst = arrset(pbm.gconst,ig_['c7'],float(-200))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xupper[ix_['Q0102']] = 1200
        pb.xupper[ix_['Q0203']] = 1200
        pb.xupper[ix_['Q0204']] = 1200
        pb.xupper[ix_['Q0305']] = 1200
        pb.xupper[ix_['Q0405']] = 1200
        pb.xupper[ix_['Q0406']] = 1200
        pb.xupper[ix_['Q0607']] = 1200
        pb.xupper[ix_['Q0705']] = 1200
        pb.xupper[ix_['Q01u0']] = 1200
        pb.xupper[ix_['y02up']] = 1200
        pb.xupper[ix_['y03up']] = 1200
        pb.xupper[ix_['y04up']] = 1200
        pb.xupper[ix_['y05up']] = 1200
        pb.xupper[ix_['y06up']] = 1200
        pb.xupper[ix_['y07up']] = 1200
        pb.xupper[ix_['yqu02']] = 1200
        pb.xupper[ix_['yqu03']] = 1200
        pb.xupper[ix_['yqu04']] = 1200
        pb.xupper[ix_['yqu05']] = 1200
        pb.xupper[ix_['yqu06']] = 1200
        pb.xupper[ix_['yqu07']] = 1200
        pb.xupper[ix_['Q0201']] = 1200
        pb.xupper[ix_['Q0302']] = 1200
        pb.xupper[ix_['Q0402']] = 1200
        pb.xupper[ix_['Q0503']] = 1200
        pb.xupper[ix_['Q0504']] = 1200
        pb.xupper[ix_['Q0604']] = 1200
        pb.xupper[ix_['Q0507']] = 1200
        pb.xupper[ix_['Q0706']] = 1200
        pb.xupper[ix_['yupu0']] = 1200
        pb.xupper[ix_['yu0uq']] = 1200
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gPOWER',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['obj0102']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0203']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0204']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0305']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0405']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0406']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0607']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        ig = ig_['obj0705']
        pbm.grftype = arrset(pbm.grftype,ig,'gPOWER')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "ONR2-MN-31-10"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gPOWER(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**2.852
        if nargout>1:
            g_ = 2.852*GVAR_**1.852
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 5.282*GVAR_**.852
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

