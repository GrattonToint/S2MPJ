from s2mpjlib import *
class  LAUNCH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The objective function to be minimized represents the total cost of
#    the development and launching of a 3 stages space launching vehicle.
#    Constraints are imposed on physical interrelations between the variables
#    and performance.
# 
#    The problem is highly non-convex. 
# 
#    Source:
#    B. Rush, J. Bracken and G. McCormick,
#    "A nonliner programming model for launch vehicle design and costing",
#    Operations Research, pp. 185-210, 1967.
# 
#    SIF input: P. Driscoll, Virginia Tech., April 1993,
#               corrected and simplified by Ph. L. Toint, May 1993.
# 
#    classification = "OOR2-MY-25-28"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LAUNCH'

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
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('AW1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AW1')
        [iv,ix_,_] = s2mpj_ii('IW1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IW1')
        [iv,ix_,_] = s2mpj_ii('MF1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MF1')
        [iv,ix_,_] = s2mpj_ii('TT1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TT1')
        [iv,ix_,_] = s2mpj_ii('PW1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PW1')
        [iv,ix_,_] = s2mpj_ii('ET1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ET1')
        [iv,ix_,_] = s2mpj_ii('S1L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'S1L')
        [iv,ix_,_] = s2mpj_ii('AW2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AW2')
        [iv,ix_,_] = s2mpj_ii('IW2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IW2')
        [iv,ix_,_] = s2mpj_ii('MF2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MF2')
        [iv,ix_,_] = s2mpj_ii('TT2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TT2')
        [iv,ix_,_] = s2mpj_ii('PW2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PW2')
        [iv,ix_,_] = s2mpj_ii('ET2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ET2')
        [iv,ix_,_] = s2mpj_ii('S2L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'S2L')
        [iv,ix_,_] = s2mpj_ii('AW3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AW3')
        [iv,ix_,_] = s2mpj_ii('IW3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IW3')
        [iv,ix_,_] = s2mpj_ii('MF3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MF3')
        [iv,ix_,_] = s2mpj_ii('TT3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TT3')
        [iv,ix_,_] = s2mpj_ii('PW3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PW3')
        [iv,ix_,_] = s2mpj_ii('ET3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ET3')
        [iv,ix_,_] = s2mpj_ii('S3L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'S3L')
        [iv,ix_,_] = s2mpj_ii('INW',ix_)
        pb.xnames=arrset(pb.xnames,iv,'INW')
        [iv,ix_,_] = s2mpj_ii('BT1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BT1')
        [iv,ix_,_] = s2mpj_ii('BT2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BT2')
        [iv,ix_,_] = s2mpj_ii('BT3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BT3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('STA1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ET1']
        pbm.A[ig,iv] = float(0.0002587)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+8))
        [ig,ig_,_] = s2mpj_ii('STA2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ET2']
        pbm.A[ig,iv] = float(0.0002587)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+8))
        [ig,ig_,_] = s2mpj_ii('STA3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ET3']
        pbm.A[ig,iv] = float(0.001958)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+8))
        [ig,ig_,_] = s2mpj_ii('INST',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['INW']
        pbm.A[ig,iv] = float(47.040096)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+8))
        [ig,ig_,_] = s2mpj_ii('LAUN',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(0.003)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(0.003)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(0.003)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(39215686.0))
        [ig,ig_,_] = s2mpj_ii('SGTH1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH1')
        iv = ix_['AW1']
        pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGTH3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH3')
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(0.6)+pbm.A[ig,iv]
        iv = ix_['AW2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGTH5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH5')
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(0.7)+pbm.A[ig,iv]
        iv = ix_['AW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGTH2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH2')
        iv = ix_['ET1']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['TT1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGTH4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH4')
        iv = ix_['ET2']
        pbm.A[ig,iv] = float(5.0)+pbm.A[ig,iv]
        iv = ix_['TT2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGTH6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SGTH6')
        iv = ix_['TT3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['ET3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI1A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SGSI1A')
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-12.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI1B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SGSI1B')
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-17.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI2A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SGSI2A')
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI2B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SGSI2B')
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-13.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI3A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SGSI3A')
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-7.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SGSI3B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SGSI3B')
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-10.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW1A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'TTIW1A')
        iv = ix_['TT1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-1.2)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW1B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'TTIW1B')
        iv = ix_['TT1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-1.4)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW2A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'TTIW2A')
        iv = ix_['TT2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-0.6)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-0.6)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-0.6)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-0.6)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-0.6)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW2B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'TTIW2B')
        iv = ix_['TT2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-0.75)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW3A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'TTIW3A')
        iv = ix_['TT3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-0.7)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-0.7)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-0.7)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('TTIW3B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'TTIW3B')
        iv = ix_['TT3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-0.9)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-0.9)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-0.9)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SMF1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SMF1')
        iv = ix_['MF1']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        iv = ix_['IW1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SMF2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SMF2')
        iv = ix_['MF2']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        iv = ix_['IW2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SMF3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SMF3')
        iv = ix_['MF3']
        pbm.A[ig,iv] = float(20.0)+pbm.A[ig,iv]
        iv = ix_['IW3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['INW']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI1A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SI1A')
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(-240.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI1B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SI1B')
        iv = ix_['PW1']
        pbm.A[ig,iv] = float(-290.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI2A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SI2A')
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-240.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI2B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SI2B')
        iv = ix_['PW2']
        pbm.A[ig,iv] = float(-290.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI3A',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'SI3A')
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-340.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SI3B',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'SI3B')
        iv = ix_['PW3']
        pbm.A[ig,iv] = float(-375.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GLGCON',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'GLGCON')
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
        pbm.gconst = arrset(pbm.gconst,ig_['STA1'],float(247.963))
        pbm.gconst = arrset(pbm.gconst,ig_['STA2'],float(247.963))
        pbm.gconst = arrset(pbm.gconst,ig_['STA3'],float(32.591))
        pbm.gconst = arrset(pbm.gconst,ig_['INST'],float(35.5))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW1A'],float(24.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW1B'],float(28.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW2A'],float(12.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW2B'],float(15.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW3A'],float(14.0))
        pbm.gconst = arrset(pbm.gconst,ig_['TTIW3B'],float(18.0))
        pbm.gconst = arrset(pbm.gconst,ig_['SMF1'],float(20.0))
        pbm.gconst = arrset(pbm.gconst,ig_['SMF2'],float(20.0))
        pbm.gconst = arrset(pbm.gconst,ig_['SMF3'],float(20.0))
        pbm.gconst = arrset(pbm.gconst,ig_['GLGCON'],float(35000.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        grange = arrset(grange,ig_['GLGCON'],float(15000.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),1.0e-8)
        pb.xupper = np.full((pb.n,1),1.0e+4)
        pb.xlower[ix_['S1L']] = 125.0
        pb.xupper[ix_['S1L']] = 150.0
        pb.xlower[ix_['S2L']] = 75.0
        pb.xupper[ix_['S2L']] = 100.0
        pb.xlower[ix_['S3L']] = 50.0
        pb.xupper[ix_['S3L']] = 70.0
        pb.xlower[ix_['MF1']] = 0.25
        pb.xupper[ix_['MF1']] = 0.30
        pb.xlower[ix_['MF2']] = 0.24
        pb.xupper[ix_['MF2']] = 0.29
        pb.xlower[ix_['MF3']] = 0.16
        pb.xupper[ix_['MF3']] = 0.21
        pb.xlower[ix_['INW']] = 2.5
        pb.xupper[ix_['INW']] = 4.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['AW1']] = float(68.0)
        pb.x0[ix_['IW1']] = float(136.0)
        pb.x0[ix_['MF1']] = float(0.29988744)
        pb.x0[ix_['TT1']] = float(3733.0)
        pb.x0[ix_['PW1']] = float(2177.0)
        pb.x0[ix_['ET1']] = float(746.6)
        pb.x0[ix_['S1L']] = float(125.0)
        pb.x0[ix_['AW2']] = float(28.2)
        pb.x0[ix_['IW2']] = float(47.0)
        pb.x0[ix_['MF2']] = float(0.28939109)
        pb.x0[ix_['TT2']] = float(480.0)
        pb.x0[ix_['PW2']] = float(566.0)
        pb.x0[ix_['ET2']] = float(96.0)
        pb.x0[ix_['S2L']] = float(75.0)
        pb.x0[ix_['AW3']] = float(11.2)
        pb.x0[ix_['IW3']] = float(16.0)
        pb.x0[ix_['MF3']] = float(0.20980926)
        pb.x0[ix_['TT3']] = float(129.0)
        pb.x0[ix_['PW3']] = float(145.0)
        pb.x0[ix_['ET3']] = float(129.0)
        pb.x0[ix_['S3L']] = float(50.0)
        pb.x0[ix_['INW']] = float(2.5)
        pb.x0[ix_['BT1']] = float(155.0)
        pb.x0[ix_['BT2']] = float(314.0)
        pb.x0[ix_['BT3']] = float(403.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD1', iet_)
        elftv = loaset(elftv,it,0,'VA')
        elftv = loaset(elftv,it,1,'VB')
        elftv = loaset(elftv,it,2,'VC')
        elftv = loaset(elftv,it,3,'VD')
        elftv = loaset(elftv,it,4,'VE')
        [it,iet_,_] = s2mpj_ii( 'ePOWER', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftp = []
        elftp = loaset(elftp,it,0,'PWR')
        elftp = loaset(elftp,it,1,'SC')
        [it,iet_,_] = s2mpj_ii( 'ePROD2', iet_)
        elftv = loaset(elftv,it,0,'VA')
        elftv = loaset(elftv,it,1,'VB')
        elftv = loaset(elftv,it,2,'VC')
        elftv = loaset(elftv,it,3,'VD')
        [it,iet_,_] = s2mpj_ii( 'eX7Y', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y1')
        elftv = loaset(elftv,it,2,'Y2')
        elftv = loaset(elftv,it,3,'Y3')
        elftv = loaset(elftv,it,4,'Y4')
        elftv = loaset(elftv,it,5,'Y5')
        elftv = loaset(elftv,it,6,'Y6')
        elftv = loaset(elftv,it,7,'Y7')
        [it,iet_,_] = s2mpj_ii( 'eX5Y', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y1')
        elftv = loaset(elftv,it,2,'Y2')
        elftv = loaset(elftv,it,3,'Y3')
        elftv = loaset(elftv,it,4,'Y4')
        elftv = loaset(elftv,it,5,'Y5')
        [it,iet_,_] = s2mpj_ii( 'eX3Y', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y1')
        elftv = loaset(elftv,it,2,'Y2')
        elftv = loaset(elftv,it,3,'Y3')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eBIG1', iet_)
        elftv = loaset(elftv,it,0,'LH')
        elftv = loaset(elftv,it,1,'TH')
        elftv = loaset(elftv,it,2,'LL')
        elftv = loaset(elftv,it,3,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'XPROD1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD1')
        ielftype = arrset(ielftype, ie, iet_["ePROD1"])
        vname = 'AW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TT1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VE')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XPF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.146))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPG'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.648))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPROD2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype, ie, iet_["ePROD2"])
        vname = 'AW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'S1L'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XPL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.736))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.229))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPROD3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD1')
        ielftype = arrset(ielftype, ie, iet_["ePROD1"])
        vname = 'AW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TT2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VE')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'X2PF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.146))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'X2PG'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.648))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPROD4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype, ie, iet_["ePROD2"])
        vname = 'AW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'S2L'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'X2PL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.736))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'X2PM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-0.229))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPROD5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD1')
        ielftype = arrset(ielftype, ie, iet_["ePROD1"])
        vname = 'AW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TT3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VE')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XQF'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.539))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XQG'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.772))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1000.0))
        ename = 'XPROD6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype, ie, iet_["ePROD2"])
        vname = 'AW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VC')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'S3L'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='VD')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XQL'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.33))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(100.0))
        ename = 'XQM'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePOWER')
        ielftype = arrset(ielftype, ie, iet_["ePOWER"])
        vname = 'ET3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='PWR')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.498))
        posep = find(elftp[ielftype[ie]],lambda x:x=='SC')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(100.0))
        ename = 'SMFE1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eX7Y')
        ielftype = arrset(ielftype, ie, iet_["eX7Y"])
        vname = 'MF1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'INW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'SMFE2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eX5Y')
        ielftype = arrset(ielftype, ie, iet_["eX5Y"])
        vname = 'MF2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'INW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'SMFE3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eX3Y')
        ielftype = arrset(ielftype, ie, iet_["eX3Y"])
        vname = 'MF3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'IW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'INW'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'TT1BT1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'TT1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'TT2BT2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'TT2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'TT3BT3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'TT3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XBIG11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBIG1')
        ielftype = arrset(ielftype, ie, iet_["eBIG1"])
        vname = 'TT1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LL')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XBIG12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBIG1')
        ielftype = arrset(ielftype, ie, iet_["eBIG1"])
        vname = 'TT2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LL')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'XBIG13'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eBIG1')
        ielftype = arrset(ielftype, ie, iet_["eBIG1"])
        vname = 'TT3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'BT3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TH')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PW3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='LL')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MF3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,1.0e-8,1.0e+4,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSUMM',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['STA1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5272.77))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(160.909))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPG'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(282.874))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.64570846))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPL'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(31.136196))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.092112))
        ig = ig_['STA2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5272.77))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X2PF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(160.909))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X2PG'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(282.874))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.64570846))
        ig = ig_['STA1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X2PL'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(31.136196))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X2PM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.092112))
        ig = ig_['STA3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(5272.77))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XQF'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(181.806))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XQG'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(232.57))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XPROD6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.49783215))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XQL'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.22424514))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XQM'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(20.708238))
        ig = ig_['LAUN']
        pbm.grftype = arrset(pbm.grftype,ig,'gSUMM')
        ig = ig_['SMF1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SMFE1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SMF2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SMFE2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SMF3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['SMFE3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI1A']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT1BT1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI1B']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT1BT1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI2A']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT2BT2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI2B']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT2BT2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI3A']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT3BT3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['SI3B']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['TT3BT3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['GLGCON']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XBIG11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-32.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XBIG12'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-32.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XBIG13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-32.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle)] = grange[legrps]
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        pb.cupper[np.arange(pb.nge)] = grange[gegrps]
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-MY-25-28"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eBIG1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LG = np.log(EV_[3])
        f_   = (EV_[0]*EV_[1]*LG)/EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[3] = (EV_[0]*EV_[1])/(EV_[3]*EV_[2])
            g_[0] = (EV_[1]*LG)/EV_[2]
            g_[1] = (EV_[0]*LG)/EV_[2]
            g_[2] = -(EV_[0]*EV_[1]*LG)/(EV_[2]**2)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[3,3] = -(EV_[0]*EV_[1])/(EV_[2]*EV_[3]**2)
                H_[3,0] = EV_[1]/(EV_[3]*EV_[2])
                H_[0,3] = H_[3,0]
                H_[3,1] = EV_[0]/(EV_[3]*EV_[2])
                H_[1,3] = H_[3,1]
                H_[3,2] = -(EV_[0]*EV_[1])/(EV_[2]**2*EV_[3])
                H_[2,3] = H_[3,2]
                H_[0,1] = LG/EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = -EV_[1]*LG/EV_[2]**2
                H_[2,0] = H_[0,2]
                H_[1,2] = -EV_[0]*LG/EV_[2]**2
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*(EV_[0]*EV_[1]*LG)/(EV_[2]**3.0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EA = 1.2781
        VA0 = EV_[0]**EA
        VA1 = EA*EV_[0]**(EA-1.0)
        VA2 = EA*(EA-1.0)*EV_[0]**(EA-2.0)
        EB = -0.1959
        VB0 = EV_[1]**EB
        VB1 = EB*EV_[1]**(EB-1.0)
        VB2 = EB*(EB-1.0)*EV_[1]**(EB-2.0)
        EC = 2.4242
        VC0 = EV_[2]**EC
        VC1 = EC*EV_[2]**(EC-1.0)
        VC2 = EC*(EC-1.0)*EV_[2]**(EC-2.0)
        ED = 0.38745
        VD0 = EV_[3]**ED
        VD1 = ED*EV_[3]**(ED-1.0)
        VD2 = ED*(ED-1.0)*EV_[3]**(ED-2.0)
        EE = 0.9904
        VE0 = EV_[4]**EE
        VE1 = EE*EV_[4]**(EE-1.0)
        VE2 = EE*(EE-1.0)*EV_[4]**(EE-2.0)
        f_   = VA0*VB0*VC0*VD0*VE0
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = VA1*VB0*VC0*VD0*VE0
            g_[1] = VA0*VB1*VC0*VD0*VE0
            g_[2] = VA0*VB0*VC1*VD0*VE0
            g_[3] = VA0*VB0*VC0*VD1*VE0
            g_[4] = VA0*VB0*VC0*VD0*VE1
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = VA2*VB0*VC0*VD0*VE0
                H_[0,1] = VA1*VB1*VC0*VD0*VE0
                H_[1,0] = H_[0,1]
                H_[0,2] = VA1*VB0*VC1*VD0*VE0
                H_[2,0] = H_[0,2]
                H_[0,3] = VA1*VB0*VC0*VD1*VE0
                H_[3,0] = H_[0,3]
                H_[0,4] = VA1*VB0*VC0*VD0*VE1
                H_[4,0] = H_[0,4]
                H_[1,1] = VA0*VB2*VC0*VD0*VE0
                H_[1,2] = VA0*VB1*VC1*VD0*VE0
                H_[2,1] = H_[1,2]
                H_[1,3] = VA0*VB1*VC0*VD1*VE0
                H_[3,1] = H_[1,3]
                H_[1,4] = VA0*VB1*VC0*VD0*VE1
                H_[4,1] = H_[1,4]
                H_[2,2] = VA0*VB0*VC2*VD0*VE0
                H_[2,3] = VA0*VB0*VC1*VD1*VE0
                H_[3,2] = H_[2,3]
                H_[2,4] = VA0*VB0*VC1*VD0*VE1
                H_[4,2] = H_[2,4]
                H_[3,3] = VA0*VB0*VC0*VD2*VE0
                H_[3,4] = VA0*VB0*VC0*VD1*VE1
                H_[4,3] = H_[3,4]
                H_[4,4] = VA0*VB0*VC0*VD0*VE2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EA = 0.3322
        VA0 = EV_[0]**EA
        VA1 = EA*EV_[0]**(EA-1.0)
        VA2 = EA*(EA-1.0)*EV_[0]**(EA-2.0)
        EB = -1.5935
        VB0 = EV_[1]**EB
        VB1 = EB*EV_[1]**(EB-1.0)
        VB2 = EB*(EB-1.0)*EV_[1]**(EB-2.0)
        EC = 0.2363
        VC0 = EV_[2]**EC
        VC1 = EC*EV_[2]**(EC-1.0)
        VC2 = EC*(EC-1.0)*EV_[2]**(EC-2.0)
        ED = 0.1079
        VD0 = EV_[3]**ED
        VD1 = ED*EV_[3]**(ED-1.0)
        VD2 = ED*(ED-1.0)*EV_[3]**(ED-2.0)
        f_   = VA0*VB0*VC0*VD0
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = VA1*VB0*VC0*VD0
            g_[1] = VA0*VB1*VC0*VD0
            g_[2] = VA0*VB0*VC1*VD0
            g_[3] = VA0*VB0*VC0*VD1
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = VA2*VB0*VC0*VD0
                H_[0,1] = VA1*VB1*VC0*VD0
                H_[1,0] = H_[0,1]
                H_[0,2] = VA1*VB0*VC1*VD0
                H_[2,0] = H_[0,2]
                H_[0,3] = VA1*VB0*VC0*VD1
                H_[3,0] = H_[0,3]
                H_[1,1] = VA0*VB2*VC0*VD0
                H_[1,2] = VA0*VB1*VC1*VD0
                H_[2,1] = H_[1,2]
                H_[1,3] = VA0*VB1*VC0*VD1
                H_[3,1] = H_[1,3]
                H_[2,2] = VA0*VB0*VC2*VD0
                H_[2,3] = VA0*VB0*VC1*VD1
                H_[3,2] = H_[2,3]
                H_[3,3] = VA0*VB0*VC0*VD2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOWER(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SCPWR = pbm.elpar[iel_][0]/(pbm.elpar[iel_][1]**pbm.elpar[iel_][0])
        f_   = (EV_[0]/pbm.elpar[iel_][1])**pbm.elpar[iel_][0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SCPWR*EV_[0]**(pbm.elpar[iel_][0]-1.0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = SCPWR*(pbm.elpar[iel_][0]-1.0)*EV_[0]**(pbm.elpar[iel_][0]-2.0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eX7Y(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,8))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]+1
        U_[1,6] = U_[1,6]+1
        U_[1,7] = U_[1,7]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eX5Y(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,6))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eX3Y(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2PR(pbm,nargout,*args):

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSUMM(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**0.460
        if nargout>1:
            g_ = 0.460*GVAR_**(-0.540)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.2484*GVAR_**(-1.540)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

