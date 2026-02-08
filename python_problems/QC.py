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
#    classification = "C-COLR2-MY-9-4"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'QC'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('F1',ix_)
        self.xnames=arrset(self.xnames,iv,'F1')
        [iv,ix_,_] = s2mpj_ii('F2',ix_)
        self.xnames=arrset(self.xnames,iv,'F2')
        [iv,ix_,_] = s2mpj_ii('PBGBG',ix_)
        self.xnames=arrset(self.xnames,iv,'PBGBG')
        [iv,ix_,_] = s2mpj_ii('PBUBG',ix_)
        self.xnames=arrset(self.xnames,iv,'PBUBG')
        [iv,ix_,_] = s2mpj_ii('PGBGB',ix_)
        self.xnames=arrset(self.xnames,iv,'PGBGB')
        [iv,ix_,_] = s2mpj_ii('PUBGB',ix_)
        self.xnames=arrset(self.xnames,iv,'PUBGB')
        [iv,ix_,_] = s2mpj_ii('PBBBB',ix_)
        self.xnames=arrset(self.xnames,iv,'PBBBB')
        [iv,ix_,_] = s2mpj_ii('PUBBB',ix_)
        self.xnames=arrset(self.xnames,iv,'PUBBB')
        [iv,ix_,_] = s2mpj_ii('PBUBB',ix_)
        self.xnames=arrset(self.xnames,iv,'PBUBB')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F1']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['2'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F1']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['3'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F2']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['4'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F2']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['5'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBG']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['6'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBGB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['7'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F2']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['8'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PGBGB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['9'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBGBG']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['10'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBBBB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['11'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PGBGB']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBGB']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['12'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBGBG']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBG']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['13'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBBBB']])
        valA = np.append(valA,float(-1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBBB']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['13'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBB']])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['14'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBG']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['15'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['16'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBGB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('OBJ'+str(int(v_['17'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBBB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CON0',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON0')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBGBG']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBG']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PGBGB']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBGB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'CON2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBBBB']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PUBBB']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['PBUBB']])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F1']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['F2']])
        valA = np.append(valA,float(1.0))
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
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['1']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['3']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['7']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['11']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['12']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['OBJ'+str(int(v_['13']))],float(-1.0))
        self.gconst = arrset(self.gconst,ig_['CON0'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CON1'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CON2'],float(1.0))
        self.gconst = arrset(self.gconst,ig_['CON3'],float(v_['RMM']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xupper = np.full((self.n,1),1.0)
        self.xlower = np.zeros((self.n,1))
        self.xlower[ix_['F1']] = v_['LF1']
        self.xupper[ix_['F1']] = v_['UF1']
        self.xlower[ix_['F2']] = v_['LF2']
        self.xupper[ix_['F2']] = v_['UF2']
        self.xlower[ix_['PBGBG']] = v_['LBGBG']
        self.xupper[ix_['PBGBG']] = v_['UBGBG']
        self.xlower[ix_['PBUBG']] = v_['LBUBG']
        self.xupper[ix_['PBUBG']] = v_['UBUBG']
        self.xlower[ix_['PGBGB']] = v_['LGBGB']
        self.xupper[ix_['PGBGB']] = v_['UGBGB']
        self.xlower[ix_['PUBGB']] = v_['LUBGB']
        self.xupper[ix_['PUBGB']] = v_['UUBGB']
        self.xlower[ix_['PBBBB']] = v_['LBBBB']
        self.xupper[ix_['PBBBB']] = v_['UBBBB']
        self.xlower[ix_['PBUBB']] = 0.0
        self.xupper[ix_['PBUBB']] = 0.0
        self.xlower[ix_['PUBBB']] = 0.0
        self.xupper[ix_['PUBBB']] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.5))
        self.x0[ix_['F1']] = float(v_['SF1'])
        self.x0[ix_['F2']] = float(v_['SF2'])
        self.x0[ix_['PBGBG']] = float(v_['SBGBG'])
        self.x0[ix_['PBUBG']] = float(v_['SBUBG'])
        self.x0[ix_['PGBGB']] = float(v_['SGBGB'])
        self.x0[ix_['PUBGB']] = float(v_['SUBGB'])
        self.x0[ix_['PBBBB']] = float(v_['SBBBB'])
        self.x0[ix_['PBUBB']] = float(v_['ZERO'])
        self.x0[ix_['PUBBB']] = float(v_['ZERO'])
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype,ie,iet_["en2PROD"])
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype,ie,iet_["en2PROD"])
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype,ie,iet_["en2PROD"])
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['3']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en2PROD')
        ielftype = arrset(ielftype,ie,iet_["en2PROD"])
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['4']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBBB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eI2PROD')
        ielftype = arrset(ielftype,ie,iet_["eI2PROD"])
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBGBG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['5']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eI3PROD')
        ielftype = arrset(ielftype,ie,iet_["eI3PROD"])
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBGBG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['6']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBG'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eI2PROD')
        ielftype = arrset(ielftype,ie,iet_["eI2PROD"])
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PGBGB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['7']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eI3PROD')
        ielftype = arrset(ielftype,ie,iet_["eI3PROD"])
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PGBGB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['8']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBGB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'en3PRODI')
        ielftype = arrset(ielftype,ie,iet_["en3PRODI"])
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'F2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBBBB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PUBBB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E'+str(int(v_['9']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'PBUBB'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,float(1.0),float(1.5))
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        [it,igt_,_] = s2mpj_ii('gLOG',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        ig = ig_['OBJ'+str(int(v_['1']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RS1']))
        ig = ig_['OBJ'+str(int(v_['2']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RS2']))
        ig = ig_['OBJ'+str(int(v_['3']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RS3']))
        ig = ig_['OBJ'+str(int(v_['4']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RS4']))
        ig = ig_['OBJ'+str(int(v_['5']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBU']))
        ig = ig_['OBJ'+str(int(v_['6']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNUB']))
        ig = ig_['OBJ'+str(int(v_['7']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['5']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['6']))])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['7']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['8']))])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['9']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNUU']))
        ig = ig_['OBJ'+str(int(v_['8']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNGB']))
        ig = ig_['OBJ'+str(int(v_['9']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBG']))
        ig = ig_['OBJ'+str(int(v_['10']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBB']))
        ig = ig_['OBJ'+str(int(v_['11']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNGBUU']))
        ig = ig_['OBJ'+str(int(v_['12']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBGUU']))
        ig = ig_['OBJ'+str(int(v_['13']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBBUU']))
        ig = ig_['OBJ'+str(int(v_['14']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBGBU']))
        ig = ig_['OBJ'+str(int(v_['15']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBBBU']))
        ig = ig_['OBJ'+str(int(v_['16']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNGBUB']))
        ig = ig_['OBJ'+str(int(v_['17']))]
        self.grftype = arrset(self.grftype,ig,'gLOG')
        posgp = np.where(grftp[igt_[self.grftype[ig]]]=='P')[0]
        self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['RNBBUB']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1138.416240
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COLR2-MY-9-4"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
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
    def eI2PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]-1
        U_[1,2] = U_[1,2]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        f_   = IV_[0]*IV_[1]
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
    def eI3PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*(1.0+IV_[2])
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
    def en3PRODI(self, nargout,*args):

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
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        f_   = IV_[0]*IV_[1]*(1.0+IV_[2])
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
    def gLOG(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        T = np.absolute(GVAR_)
        SMALL = 1.0e-10
        LARGE = 1.0e+10
        ARG0 = T<=SMALL
        if ARG0!=0:
            FF = self.grpar[igr_][0]*np.log(SMALL)
        if ARG0==0:
            FF = self.grpar[igr_][0]*np.log(T)
        if ARG0!=0:
            GG = self.grpar[igr_][0]*LARGE
        if ARG0==0:
            GG = self.grpar[igr_][0]/T
        if ARG0!=0:
            HH = -self.grpar[igr_][0]*LARGE**2
        if ARG0==0:
            HH = -self.grpar[igr_][0]/T**2
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

