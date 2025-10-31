from s2mpjlib import *
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
#    classification = "C-CONR2-MN-31-10"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'WATER'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [ig,ig_,_] = s2mpj_ii('obj0102',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(105665.6))
        [ig,ig_,_] = s2mpj_ii('obj0203',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(3613.412))
        [ig,ig_,_] = s2mpj_ii('obj0204',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(105665.6))
        [ig,ig_,_] = s2mpj_ii('obj0305',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(890.1553))
        [ig,ig_,_] = s2mpj_ii('obj0405',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(76.66088))
        [ig,ig_,_] = s2mpj_ii('obj0406',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(55145.82))
        [ig,ig_,_] = s2mpj_ii('obj0607',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(26030.46))
        [ig,ig_,_] = s2mpj_ii('obj0705',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(890.1553))
        [ig,ig_,_] = s2mpj_ii('obj',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('c1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c1')
        [ig,ig_,_] = s2mpj_ii('c2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c2')
        [ig,ig_,_] = s2mpj_ii('c3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c3')
        [ig,ig_,_] = s2mpj_ii('c4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c4')
        [ig,ig_,_] = s2mpj_ii('c5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c5')
        [ig,ig_,_] = s2mpj_ii('c6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c6')
        [ig,ig_,_] = s2mpj_ii('c7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c7')
        [ig,ig_,_] = s2mpj_ii('c8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c8')
        [ig,ig_,_] = s2mpj_ii('c9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c9')
        [ig,ig_,_] = s2mpj_ii('c10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'c10')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('Q0102',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0102')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0102']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c1']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0102',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0102')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0203',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0203')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0203']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0203',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0203')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0204',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0204')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0204']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0204',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0204')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0305',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0305')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0305']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0305',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0305')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0405',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0405')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0405']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0405',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0405')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0406',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0406')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0406']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0406',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0406')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0607',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0607')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0607']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0607',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0607')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0705',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0705')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj0705']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0705',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0705')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q01u0',ix_)
        self.xnames=arrset(self.xnames,iv,'Q01u0')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c1']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q01u0',ix_)
        self.xnames=arrset(self.xnames,iv,'Q01u0')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c8']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y02up',ix_)
        self.xnames=arrset(self.xnames,iv,'y02up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y02up',ix_)
        self.xnames=arrset(self.xnames,iv,'y02up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y03up',ix_)
        self.xnames=arrset(self.xnames,iv,'y03up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y03up',ix_)
        self.xnames=arrset(self.xnames,iv,'y03up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y04up',ix_)
        self.xnames=arrset(self.xnames,iv,'y04up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y04up',ix_)
        self.xnames=arrset(self.xnames,iv,'y04up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y05up',ix_)
        self.xnames=arrset(self.xnames,iv,'y05up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y05up',ix_)
        self.xnames=arrset(self.xnames,iv,'y05up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y06up',ix_)
        self.xnames=arrset(self.xnames,iv,'y06up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y06up',ix_)
        self.xnames=arrset(self.xnames,iv,'y06up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('y07up',ix_)
        self.xnames=arrset(self.xnames,iv,'y07up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(210))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('y07up',ix_)
        self.xnames=arrset(self.xnames,iv,'y07up')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu02',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu02')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-175))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu02',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu02')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yqu03',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu03')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-190))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu03',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu03')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yqu04',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu04')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-185))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu04',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu04')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yqu05',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu05')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-180))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu05',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu05')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yqu06',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu06')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-195))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu06',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu06')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yqu07',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu07')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['obj']])
        valA = np.append(valA,float(-190))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('yqu07',ix_)
        self.xnames=arrset(self.xnames,iv,'yqu07')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0201',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0201')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c1']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0302',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0302')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0402',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0402')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c2']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0503',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0503')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c3']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0504',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0504')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0604',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0604')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c4']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('Q0507',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0507')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c5']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(-1))
        [iv,ix_,_] = s2mpj_ii('Q0706',ix_)
        self.xnames=arrset(self.xnames,iv,'Q0706')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c6']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c7']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yupu0',ix_)
        self.xnames=arrset(self.xnames,iv,'yupu0')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c8']])
        valA = np.append(valA,float(-1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c9']])
        valA = np.append(valA,float(1))
        [iv,ix_,_] = s2mpj_ii('yu0uq',ix_)
        self.xnames=arrset(self.xnames,iv,'yu0uq')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c8']])
        valA = np.append(valA,float(1))
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['c10']])
        valA = np.append(valA,float(-1))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
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
        self.gconst = arrset(self.gconst,ig_['c1'],float(1120))
        self.gconst = arrset(self.gconst,ig_['c2'],float(-100))
        self.gconst = arrset(self.gconst,ig_['c3'],float(-100))
        self.gconst = arrset(self.gconst,ig_['c4'],float(-120))
        self.gconst = arrset(self.gconst,ig_['c5'],float(-270))
        self.gconst = arrset(self.gconst,ig_['c6'],float(-330))
        self.gconst = arrset(self.gconst,ig_['c7'],float(-200))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['Q0102']] = 1200
        self.xupper[ix_['Q0203']] = 1200
        self.xupper[ix_['Q0204']] = 1200
        self.xupper[ix_['Q0305']] = 1200
        self.xupper[ix_['Q0405']] = 1200
        self.xupper[ix_['Q0406']] = 1200
        self.xupper[ix_['Q0607']] = 1200
        self.xupper[ix_['Q0705']] = 1200
        self.xupper[ix_['Q01u0']] = 1200
        self.xupper[ix_['y02up']] = 1200
        self.xupper[ix_['y03up']] = 1200
        self.xupper[ix_['y04up']] = 1200
        self.xupper[ix_['y05up']] = 1200
        self.xupper[ix_['y06up']] = 1200
        self.xupper[ix_['y07up']] = 1200
        self.xupper[ix_['yqu02']] = 1200
        self.xupper[ix_['yqu03']] = 1200
        self.xupper[ix_['yqu04']] = 1200
        self.xupper[ix_['yqu05']] = 1200
        self.xupper[ix_['yqu06']] = 1200
        self.xupper[ix_['yqu07']] = 1200
        self.xupper[ix_['Q0201']] = 1200
        self.xupper[ix_['Q0302']] = 1200
        self.xupper[ix_['Q0402']] = 1200
        self.xupper[ix_['Q0503']] = 1200
        self.xupper[ix_['Q0504']] = 1200
        self.xupper[ix_['Q0604']] = 1200
        self.xupper[ix_['Q0507']] = 1200
        self.xupper[ix_['Q0706']] = 1200
        self.xupper[ix_['yupu0']] = 1200
        self.xupper[ix_['yu0uq']] = 1200
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gPOWER',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['obj0102']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0203']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0204']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0305']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0405']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0406']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0607']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        ig = ig_['obj0705']
        self.grftype = arrset(self.grftype,ig,'gPOWER')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION           1.054938D+04
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CONR2-MN-31-10"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gPOWER(self,nargout,*args):

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

