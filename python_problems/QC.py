from s2mpjlib import *
class  QC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QC
#    *********
# 
#    Source: Quality Control problem 104 from
#    Betty Schultz and Ben Reiser.
# 
#    SIF input: Andrew Conn, August 1992.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "OLR2-MY-9-4"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'QC'

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
        v_['M'] = 9
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['15'] = 15
        v_['16'] = 16
        v_['17'] = 17
        v_['N'] = 10048
        v_['NGG'] = 9900
        v_['NGB'] = 35
        v_['NBG'] = 15
        v_['NBU'] = 18
        v_['NBB'] = 2
        v_['NUU'] = 40
        v_['NUB'] = 21
        v_['NGBUU'] = 9
        v_['NBGUU'] = 2
        v_['NBGBU'] = 2
        v_['NGBUB'] = 3
        v_['NBBUU'] = 1
        v_['NBBBU'] = 0
        v_['NBBUB'] = 0
        v_['MM'] = v_['N']-v_['NGG']
        v_['S1'] = v_['NGG']+v_['NGB']
        v_['S1'] = v_['S1']+v_['NGBUU']
        v_['S1'] = v_['S1']+v_['NGBUB']
        v_['S2'] = v_['NBG']+v_['NBB']
        v_['S2'] = v_['S2']+v_['NBU']
        v_['S2'] = v_['S2']+v_['NBGUU']
        v_['S2'] = v_['S2']+v_['NBBUU']
        v_['S2'] = v_['S2']+v_['NBGBU']
        v_['S2'] = v_['S2']+v_['NBBBU']
        v_['S2'] = v_['S2']+v_['NBBUB']
        v_['S3'] = v_['NGG']+v_['NBG']
        v_['S3'] = v_['S3']+v_['NBGUU']
        v_['S3'] = v_['S3']+v_['NBGBU']
        v_['S4'] = v_['NGB']+v_['NBB']
        v_['S4'] = v_['S4']+v_['NUB']
        v_['S4'] = v_['S4']+v_['NGBUU']
        v_['S4'] = v_['S4']+v_['NBBUU']
        v_['S4'] = v_['S4']+v_['NBBBU']
        v_['S4'] = v_['S4']+v_['NGBUB']
        v_['S4'] = v_['S4']+v_['NBBUB']
        v_['U1'] = v_['NGB']+v_['NGBUU']
        v_['U1'] = v_['U1']+v_['NGBUB']
        v_['L1'] = v_['U1']+v_['NUU']
        v_['L1'] = v_['L1']+v_['NUB']
        v_['U2'] = v_['NBG']+v_['NBGUU']
        v_['U2'] = v_['U2']+v_['NBGBU']
        v_['L2'] = v_['U2']+v_['NUU']
        v_['L2'] = v_['L2']+v_['NBU']
        v_['U3'] = v_['NBB']+v_['NBBUU']
        v_['U3'] = v_['U3']+v_['NBBBU']
        v_['U3'] = v_['U3']+v_['NBBUB']
        v_['L3'] = v_['U3']+v_['NUU']
        v_['L3'] = v_['L3']+v_['NBU']
        v_['L3'] = v_['L3']+v_['NUB']
        v_['U4'] = v_['U2']+v_['NBU']
        v_['L4'] = v_['U4']+v_['NUU']
        v_['U5'] = v_['U1']+v_['NUB']
        v_['L5'] = v_['U5']+v_['NUU']
        v_['U6'] = v_['U3']+v_['NUB']
        v_['L6'] = v_['U6']+v_['NBU']
        v_['L6'] = v_['L6']+v_['NUU']
        v_['U7'] = v_['U3']+v_['NBU']
        v_['L7'] = v_['U7']+v_['NUB']
        v_['L7'] = v_['L7']+v_['NUU']
        v_['L8'] = v_['U2']+v_['NBBUU']
        v_['L8'] = v_['L8']+v_['NBBBU']
        v_['L8'] = v_['L8']+v_['NBBUB']
        v_['L8'] = v_['L8']+v_['NBU']
        v_['L8'] = v_['L8']+v_['NBB']
        v_['U8'] = v_['L8']+v_['NUU']
        v_['U8'] = v_['U8']+v_['NUB']
        v_['L9'] = v_['U1']+v_['NBBUU']
        v_['L9'] = v_['L9']+v_['NBBBU']
        v_['L9'] = v_['L9']+v_['NBBUB']
        v_['L9'] = v_['L9']+v_['NUB']
        v_['L9'] = v_['L9']+v_['NBB']
        v_['U9'] = v_['L9']+v_['NUU']
        v_['U9'] = v_['U9']+v_['NBU']
        v_['ZERO'] = 0.0
        v_['TWO'] = 2.0
        v_['RS1'] = float(v_['S1'])
        v_['RS2'] = float(v_['S2'])
        v_['RS3'] = float(v_['S3'])
        v_['RS4'] = float(v_['S4'])
        v_['RNBG'] = float(v_['NBG'])
        v_['RNBGBU'] = float(v_['NBGBU'])
        v_['RNBGUU'] = float(v_['NBGUU'])
        v_['RNGB'] = float(v_['NGB'])
        v_['RNGBUB'] = float(v_['NGBUB'])
        v_['RNGBUU'] = float(v_['NGBUU'])
        v_['RNBB'] = float(v_['NBB'])
        v_['RNBBUB'] = float(v_['NBBUB'])
        v_['RNBBBU'] = float(v_['NBBBU'])
        v_['RNBBUU'] = float(v_['NBBUU'])
        v_['RNBGBU'] = float(v_['NBGBU'])
        v_['RNUU'] = float(v_['NUU'])
        v_['RNUB'] = float(v_['NUB'])
        v_['RNBU'] = float(v_['NBU'])
        v_['RN'] = float(v_['N'])
        v_['RL1'] = float(v_['L1'])
        v_['RU1'] = float(v_['U1'])
        v_['RL2'] = float(v_['L2'])
        v_['RU2'] = float(v_['U2'])
        v_['RL3'] = float(v_['L3'])
        v_['RU3'] = float(v_['U3'])
        v_['RL4'] = float(v_['L4'])
        v_['RU4'] = float(v_['U4'])
        v_['RL5'] = float(v_['L5'])
        v_['RU5'] = float(v_['U5'])
        v_['RL6'] = float(v_['L6'])
        v_['RU6'] = float(v_['U6'])
        v_['RL7'] = float(v_['L7'])
        v_['RU7'] = float(v_['U7'])
        v_['RL8'] = float(v_['L8'])
        v_['RU8'] = float(v_['U8'])
        v_['RL9'] = float(v_['L9'])
        v_['RU9'] = float(v_['U9'])
        v_['LF1'] = v_['RL8']/v_['RN']
        v_['UF1'] = v_['RU8']/v_['RN']
        v_['SF1'] = v_['LF1']+v_['UF1']
        v_['SF1'] = v_['SF1']/v_['TWO']
        v_['LF2'] = v_['RL9']/v_['RN']
        v_['UF2'] = v_['RU9']/v_['RN']
        v_['SF2'] = v_['LF2']+v_['UF2']
        v_['SF2'] = v_['SF2']/v_['TWO']
        v_['LGBGB'] = v_['RNGB']/v_['RL1']
        v_['UGBGB'] = v_['RNGB']/v_['RU1']
        v_['SGBGB'] = v_['LGBGB']+v_['UGBGB']
        v_['SGBGB'] = v_['SGBGB']/v_['TWO']
        v_['LUBGB'] = v_['RNGBUB']/v_['RL5']
        v_['UUBGB'] = v_['RNGBUB']+v_['RNUB']
        v_['UUBGB'] = v_['UUBGB']/v_['RU5']
        v_['SUBGB'] = v_['LUBGB']+v_['UUBGB']
        v_['SUBGB'] = v_['SUBGB']/v_['TWO']
        v_['LBGBG'] = v_['RNBG']/v_['RL2']
        v_['UBGBG'] = v_['RNBG']/v_['RU2']
        v_['SBGBG'] = v_['LBGBG']+v_['UBGBG']
        v_['SBGBG'] = v_['SBGBG']/v_['TWO']
        v_['LBUBG'] = v_['RNBGBU']/v_['RL4']
        v_['UBUBG'] = v_['RNBGBU']+v_['RNBU']
        v_['UBUBG'] = v_['UBUBG']/v_['RU4']
        v_['SBUBG'] = v_['LBUBG']+v_['UBUBG']
        v_['SBUBG'] = v_['SBUBG']/v_['TWO']
        v_['LBBBB'] = v_['RNBB']/v_['RL3']
        v_['UBBBB'] = v_['RNBB']/v_['RU3']
        v_['SBBBB'] = v_['LBBBB']+v_['UBBBB']
        v_['SBBBB'] = v_['SBBBB']/v_['TWO']
        v_['LBUBB'] = v_['RNBBBU']/v_['RL7']
        v_['UBUBB'] = v_['RNBBBU']+v_['RNBU']
        v_['UBUBB'] = v_['UBUBB']/v_['RU7']
        v_['SBUBB'] = v_['LBUBB']+v_['UBUBB']
        v_['SBUBB'] = v_['SBUBB']/v_['TWO']
        v_['LUBBB'] = v_['RNBBUB']/v_['RL6']
        v_['UUBBB'] = v_['RNBBUB']+v_['RNUB']
        v_['UUBBB'] = v_['UUBBB']/v_['RU6']
        v_['SUBBB'] = v_['LUBBB']+v_['UUBBB']
        v_['SUBBB'] = v_['SUBBB']/v_['TWO']
        v_['RMM'] = float(v_['MM'])
        v_['RMM'] = v_['RMM']/v_['RN']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('F1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'F1')
        [iv,ix_,_] = s2mpj_ii('F2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'F2')
        [iv,ix_,_] = s2mpj_ii('PBGBG',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PBGBG')
        [iv,ix_,_] = s2mpj_ii('PBUBG',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PBUBG')
        [iv,ix_,_] = s2mpj_ii('PGBGB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PGBGB')
        [iv,ix_,_] = s2mpj_ii('PUBGB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PUBGB')
        [iv,ix_,_] = s2mpj_ii('PBBBB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PBBBB')
        [iv,ix_,_] = s2mpj_ii('PUBBB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PUBBB')
        [iv,ix_,_] = s2mpj_ii('PBUBB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PBUBB')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['2'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['3'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['4'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['5'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBUBG']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['6'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PUBGB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['7'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['F1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['F2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['8'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PGBGB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['9'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBGBG']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['10'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBBBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['11'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PGBGB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PUBGB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['12'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBGBG']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PBUBG']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['13'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBBBB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['PUBBB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['13'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBUBB']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['14'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBUBG']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['15'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PBUBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['16'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PUBGB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['17'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PUBBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON0',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON0')
        iv = ix_['PBGBG']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['PBUBG']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON1')
        iv = ix_['PGBGB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['PUBGB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON2')
        iv = ix_['PBBBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['PUBBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['PBUBB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON3')
        iv = ix_['F1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['F2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['1']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['3']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['7']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['11']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['12']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(int(v_['13']))],float(-1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON0'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON3'],float(v_['RMM']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),1.0)
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['F1']] = v_['LF1']
        pb.xupper[ix_['F1']] = v_['UF1']
        pb.xlower[ix_['F2']] = v_['LF2']
        pb.xupper[ix_['F2']] = v_['UF2']
        pb.xlower[ix_['PBGBG']] = v_['LBGBG']
        pb.xupper[ix_['PBGBG']] = v_['UBGBG']
        pb.xlower[ix_['PBUBG']] = v_['LBUBG']
        pb.xupper[ix_['PBUBG']] = v_['UBUBG']
        pb.xlower[ix_['PGBGB']] = v_['LGBGB']
        pb.xupper[ix_['PGBGB']] = v_['UGBGB']
        pb.xlower[ix_['PUBGB']] = v_['LUBGB']
        pb.xupper[ix_['PUBGB']] = v_['UUBGB']
        pb.xlower[ix_['PBBBB']] = v_['LBBBB']
        pb.xupper[ix_['PBBBB']] = v_['UBBBB']
        pb.xlower[ix_['PBUBB']] = 0.0
        pb.xupper[ix_['PBUBB']] = 0.0
        pb.xlower[ix_['PUBBB']] = 0.0
        pb.xupper[ix_['PUBBB']] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.5))
        pb.x0[ix_['F1']] = float(v_['SF1'])
        pb.x0[ix_['F2']] = float(v_['SF2'])
        pb.x0[ix_['PBGBG']] = float(v_['SBGBG'])
        pb.x0[ix_['PBUBG']] = float(v_['SBUBG'])
        pb.x0[ix_['PGBGB']] = float(v_['SGBGB'])
        pb.x0[ix_['PUBGB']] = float(v_['SUBGB'])
        pb.x0[ix_['PBBBB']] = float(v_['SBBBB'])
        pb.x0[ix_['PBUBB']] = float(v_['ZERO'])
        pb.x0[ix_['PUBBB']] = float(v_['ZERO'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'eI2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        [it,iet_,_] = s2mpj_ii( 'eI3PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        [it,iet_,_] = s2mpj_ii( 'en3PRODI', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype, ie, iet_["en2PROD"])
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype, ie, iet_["en2PROD"])
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype, ie, iet_["en2PROD"])
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype, ie, iet_["en2PROD"])
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBBB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eI2PROD')
        ielftype = arrset(ielftype, ie, iet_["eI2PROD"])
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBGBG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eI3PROD')
        ielftype = arrset(ielftype, ie, iet_["eI3PROD"])
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBGBG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eI2PROD')
        ielftype = arrset(ielftype, ie, iet_["eI2PROD"])
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PGBGB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eI3PROD')
        ielftype = arrset(ielftype, ie, iet_["eI3PROD"])
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PGBGB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en3PRODI')
        ielftype = arrset(ielftype, ie, iet_["en3PRODI"])
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBBBB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBBB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBB'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.0,1.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        ig = ig_['OBJ'+str(int(v_['1']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RS1']))
        ig = ig_['OBJ'+str(int(v_['2']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RS2']))
        ig = ig_['OBJ'+str(int(v_['3']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RS3']))
        ig = ig_['OBJ'+str(int(v_['4']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RS4']))
        ig = ig_['OBJ'+str(int(v_['5']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBU']))
        ig = ig_['OBJ'+str(int(v_['6']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNUB']))
        ig = ig_['OBJ'+str(int(v_['7']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['5']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['6']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['7']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['8']))])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['9']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNUU']))
        ig = ig_['OBJ'+str(int(v_['8']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNGB']))
        ig = ig_['OBJ'+str(int(v_['9']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBG']))
        ig = ig_['OBJ'+str(int(v_['10']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBB']))
        ig = ig_['OBJ'+str(int(v_['11']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNGBUU']))
        ig = ig_['OBJ'+str(int(v_['12']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBGUU']))
        ig = ig_['OBJ'+str(int(v_['13']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBBUU']))
        ig = ig_['OBJ'+str(int(v_['14']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBGBU']))
        ig = ig_['OBJ'+str(int(v_['15']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBBBU']))
        ig = ig_['OBJ'+str(int(v_['16']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNGBUB']))
        ig = ig_['OBJ'+str(int(v_['17']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLOG')
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RNBBUB']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
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
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-MY-9-4"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(pbm,nargout,*args):

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
                H_[0,0] = 0.0
                H_[1,1] = 0.0
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eI2PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]-1
        U_[1,2] = U_[1,2]-1
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
                H_[0,0] = 0.0
                H_[1,1] = 0.0
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
    def eI3PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        f_   = IV_[0]*IV_[1]*(1.0+IV_[2])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*(1.0+IV_[2])
            g_[1] = IV_[0]*(1.0+IV_[2])
            g_[2] = IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 0.0
                H_[0,1] = 1.0+IV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = 0.0
                H_[1,2] = IV_[0]
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3PRODI(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,5))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]-1
        U_[2,4] = U_[2,4]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        f_   = IV_[0]*IV_[1]*(1.0+IV_[2])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*(1.0+IV_[2])
            g_[1] = IV_[0]*(1.0+IV_[2])
            g_[2] = IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 0.0
                H_[0,1] = 1.0+IV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]
                H_[2,0] = H_[0,2]
                H_[1,1] = 0.0
                H_[1,2] = IV_[0]
                H_[2,1] = H_[1,2]
                H_[2,2] = 0.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLOG(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        T = np.absolute(GVAR_)
        SMALL = 1.0e-10
        LARGE = 1.0e+10
        ARG0 = T<=SMALL
        if ARG0!=0:
            FF = pbm.grpar[igr_][0]*np.log(SMALL)
        if ARG0==0:
            FF = pbm.grpar[igr_][0]*np.log(T)
        if ARG0!=0:
            GG = pbm.grpar[igr_][0]*LARGE
        if ARG0==0:
            GG = pbm.grpar[igr_][0]/T
        if ARG0!=0:
            HH = -pbm.grpar[igr_][0]*LARGE**2
        if ARG0==0:
            HH = -pbm.grpar[igr_][0]/T**2
        f_= FF
        if nargout>1:
            g_ = GG
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = HH
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

