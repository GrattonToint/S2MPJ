from s2mpjlib import *
class  QPCBLEND(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QPCBLEND
#    *********
# 
#    Source: a variant on the BLEND linear programming problem
#    with an additional CONVEX diagonal Hessian matrix as given by
#    N. I. M. Gould, "An algorithm for large-scale quadratic programming",
#    IMA J. Num. Anal (1991), 11, 299-324, problem class 4.
# 
#    SIF input: Nick Gould, January 1993
# 
#    classification = "C-CQLR2-MN-83-74"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'QPCBLEND'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'1')
        [ig,ig_,_] = s2mpj_ii('2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'2')
        [ig,ig_,_] = s2mpj_ii('3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'3')
        [ig,ig_,_] = s2mpj_ii('4',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'4')
        [ig,ig_,_] = s2mpj_ii('5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'5')
        [ig,ig_,_] = s2mpj_ii('6',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'6')
        [ig,ig_,_] = s2mpj_ii('7',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'7')
        [ig,ig_,_] = s2mpj_ii('8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'8')
        [ig,ig_,_] = s2mpj_ii('9',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'9')
        [ig,ig_,_] = s2mpj_ii('10',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'10')
        [ig,ig_,_] = s2mpj_ii('11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'11')
        [ig,ig_,_] = s2mpj_ii('12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'12')
        [ig,ig_,_] = s2mpj_ii('13',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'13')
        [ig,ig_,_] = s2mpj_ii('14',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'14')
        [ig,ig_,_] = s2mpj_ii('15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'15')
        [ig,ig_,_] = s2mpj_ii('16',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'16')
        [ig,ig_,_] = s2mpj_ii('17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'17')
        [ig,ig_,_] = s2mpj_ii('18',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'18')
        [ig,ig_,_] = s2mpj_ii('19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'19')
        [ig,ig_,_] = s2mpj_ii('20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'20')
        [ig,ig_,_] = s2mpj_ii('21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'21')
        [ig,ig_,_] = s2mpj_ii('22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'22')
        [ig,ig_,_] = s2mpj_ii('23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'23')
        [ig,ig_,_] = s2mpj_ii('24',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'24')
        [ig,ig_,_] = s2mpj_ii('25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'25')
        [ig,ig_,_] = s2mpj_ii('26',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'26')
        [ig,ig_,_] = s2mpj_ii('27',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'27')
        [ig,ig_,_] = s2mpj_ii('28',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'28')
        [ig,ig_,_] = s2mpj_ii('29',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'29')
        [ig,ig_,_] = s2mpj_ii('30',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'30')
        [ig,ig_,_] = s2mpj_ii('31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'31')
        [ig,ig_,_] = s2mpj_ii('32',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'32')
        [ig,ig_,_] = s2mpj_ii('33',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'33')
        [ig,ig_,_] = s2mpj_ii('34',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'34')
        [ig,ig_,_] = s2mpj_ii('35',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'35')
        [ig,ig_,_] = s2mpj_ii('36',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'36')
        [ig,ig_,_] = s2mpj_ii('37',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'37')
        [ig,ig_,_] = s2mpj_ii('38',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'38')
        [ig,ig_,_] = s2mpj_ii('39',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'39')
        [ig,ig_,_] = s2mpj_ii('40',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'40')
        [ig,ig_,_] = s2mpj_ii('41',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'41')
        [ig,ig_,_] = s2mpj_ii('42',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'42')
        [ig,ig_,_] = s2mpj_ii('43',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'43')
        [ig,ig_,_] = s2mpj_ii('44',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'44')
        [ig,ig_,_] = s2mpj_ii('45',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'45')
        [ig,ig_,_] = s2mpj_ii('46',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'46')
        [ig,ig_,_] = s2mpj_ii('47',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'47')
        [ig,ig_,_] = s2mpj_ii('48',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'48')
        [ig,ig_,_] = s2mpj_ii('49',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'49')
        [ig,ig_,_] = s2mpj_ii('50',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'50')
        [ig,ig_,_] = s2mpj_ii('51',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'51')
        [ig,ig_,_] = s2mpj_ii('52',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'52')
        [ig,ig_,_] = s2mpj_ii('53',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'53')
        [ig,ig_,_] = s2mpj_ii('54',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'54')
        [ig,ig_,_] = s2mpj_ii('55',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'55')
        [ig,ig_,_] = s2mpj_ii('56',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'56')
        [ig,ig_,_] = s2mpj_ii('57',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'57')
        [ig,ig_,_] = s2mpj_ii('58',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'58')
        [ig,ig_,_] = s2mpj_ii('59',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'59')
        [ig,ig_,_] = s2mpj_ii('60',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'60')
        [ig,ig_,_] = s2mpj_ii('61',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'61')
        [ig,ig_,_] = s2mpj_ii('62',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'62')
        [ig,ig_,_] = s2mpj_ii('63',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'63')
        [ig,ig_,_] = s2mpj_ii('64',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'64')
        [ig,ig_,_] = s2mpj_ii('65',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'65')
        [ig,ig_,_] = s2mpj_ii('66',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'66')
        [ig,ig_,_] = s2mpj_ii('67',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'67')
        [ig,ig_,_] = s2mpj_ii('68',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'68')
        [ig,ig_,_] = s2mpj_ii('69',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'69')
        [ig,ig_,_] = s2mpj_ii('70',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'70')
        [ig,ig_,_] = s2mpj_ii('71',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'71')
        [ig,ig_,_] = s2mpj_ii('72',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'72')
        [ig,ig_,_] = s2mpj_ii('73',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'73')
        [ig,ig_,_] = s2mpj_ii('74',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'74')
        [ig,ig_,_] = s2mpj_ii('C',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['2']
        self.A[ig,iv] = float(-.537)+self.A[ig,iv]
        ig = ig_['3']
        self.A[ig,iv] = float(-.131)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['4']
        self.A[ig,iv] = float(-.1155)+self.A[ig,iv]
        ig = ig_['5']
        self.A[ig,iv] = float(-.0365)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['6']
        self.A[ig,iv] = float(-.143)+self.A[ig,iv]
        ig = ig_['7']
        self.A[ig,iv] = float(-.037)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['40']
        self.A[ig,iv] = float(.003)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(.0587)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['42']
        self.A[ig,iv] = float(.15)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.302)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        ig = ig_['67']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(3.2)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['1']
        self.A[ig,iv] = float(-.2931)+self.A[ig,iv]
        ig = ig_['3']
        self.A[ig,iv] = float(-.117)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['4']
        self.A[ig,iv] = float(-.0649)+self.A[ig,iv]
        ig = ig_['5']
        self.A[ig,iv] = float(-.1233)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['6']
        self.A[ig,iv] = float(-.2217)+self.A[ig,iv]
        ig = ig_['8']
        self.A[ig,iv] = float(-.18)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['39']
        self.A[ig,iv] = float(.0042)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.003)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['41']
        self.A[ig,iv] = float(.1053)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.185)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['43']
        self.A[ig,iv] = float(.384)+self.A[ig,iv]
        ig = ig_['50']
        self.A[ig,iv] = float(-.00862)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['51']
        self.A[ig,iv] = float(-.00862)+self.A[ig,iv]
        ig = ig_['56']
        self.A[ig,iv] = float(-.0101)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['57']
        self.A[ig,iv] = float(-.0101)+self.A[ig,iv]
        ig = ig_['68']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        ig = ig_['C']
        self.A[ig,iv] = float(2.87)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['9']
        self.A[ig,iv] = float(-.0277)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['10']
        self.A[ig,iv] = float(-.0563)+self.A[ig,iv]
        ig = ig_['11']
        self.A[ig,iv] = float(-.199)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['12']
        self.A[ig,iv] = float(-.6873)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.017)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['40']
        self.A[ig,iv] = float(.01303)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(.0506)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['42']
        self.A[ig,iv] = float(.209)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.495)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['9']
        self.A[ig,iv] = float(-.0112)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['10']
        self.A[ig,iv] = float(-.0378)+self.A[ig,iv]
        ig = ig_['11']
        self.A[ig,iv] = float(-.1502)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['12']
        self.A[ig,iv] = float(-.7953)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.0099)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['40']
        self.A[ig,iv] = float(.01303)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(.0448)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['42']
        self.A[ig,iv] = float(.185)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.721)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        ig = ig_['9']
        self.A[ig,iv] = float(-.175)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.27)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        ig = ig_['11']
        self.A[ig,iv] = float(-.028)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.455)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        ig = ig_['21']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.01303)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        ig = ig_['41']
        self.A[ig,iv] = float(.0506)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.209)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        ig = ig_['43']
        self.A[ig,iv] = float(.495)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        ig = ig_['9']
        self.A[ig,iv] = float(-.271)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.3285)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        ig = ig_['11']
        self.A[ig,iv] = float(-.0255)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.2656)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        ig = ig_['18']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.01303)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        ig = ig_['41']
        self.A[ig,iv] = float(.0506)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.209)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        ig = ig_['43']
        self.A[ig,iv] = float(.495)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        ig = ig_['9']
        self.A[ig,iv] = float(-.2836)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.3285)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        ig = ig_['11']
        self.A[ig,iv] = float(-.0241)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.2502)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        ig = ig_['17']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.01303)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        ig = ig_['41']
        self.A[ig,iv] = float(.0506)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.209)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        ig = ig_['43']
        self.A[ig,iv] = float(.495)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('8',ix_)
        self.xnames=arrset(self.xnames,iv,'8')
        ig = ig_['12']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['14']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('8',ix_)
        self.xnames=arrset(self.xnames,iv,'8')
        ig = ig_['39']
        self.A[ig,iv] = float(.0327)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(.094)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('8',ix_)
        self.xnames=arrset(self.xnames,iv,'8')
        ig = ig_['42']
        self.A[ig,iv] = float(.045)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.793)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('8',ix_)
        self.xnames=arrset(self.xnames,iv,'8')
        ig = ig_['C']
        self.A[ig,iv] = float(.0044)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('9',ix_)
        self.xnames=arrset(self.xnames,iv,'9')
        ig = ig_['15']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['22']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('10',ix_)
        self.xnames=arrset(self.xnames,iv,'10')
        ig = ig_['16']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['22']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('11',ix_)
        self.xnames=arrset(self.xnames,iv,'11')
        ig = ig_['14']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['15']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('12',ix_)
        self.xnames=arrset(self.xnames,iv,'12')
        ig = ig_['14']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['16']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['15']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['17']
        self.A[ig,iv] = float(-.0588)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['19']
        self.A[ig,iv] = float(-.8145)+self.A[ig,iv]
        ig = ig_['23']
        self.A[ig,iv] = float(-.0091)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['39']
        self.A[ig,iv] = float(-.8239)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.0081)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['41']
        self.A[ig,iv] = float(-.2112)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.387)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['43']
        self.A[ig,iv] = float(1.03)+self.A[ig,iv]
        ig = ig_['69']
        self.A[ig,iv] = float(1.3)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        ig = ig_['C']
        self.A[ig,iv] = float(.07)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['16']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['18']
        self.A[ig,iv] = float(-.0404)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['20']
        self.A[ig,iv] = float(-.8564)+self.A[ig,iv]
        ig = ig_['23']
        self.A[ig,iv] = float(-.0069)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['39']
        self.A[ig,iv] = float(-.7689)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.0063)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['41']
        self.A[ig,iv] = float(-.156)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.297)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['43']
        self.A[ig,iv] = float(.792)+self.A[ig,iv]
        ig = ig_['69']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        ig = ig_['C']
        self.A[ig,iv] = float(.0378)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['21']
        self.A[ig,iv] = float(-.3321)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['22']
        self.A[ig,iv] = float(-.5875)+self.A[ig,iv]
        ig = ig_['23']
        self.A[ig,iv] = float(-.362)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['39']
        self.A[ig,iv] = float(2.3)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.2049)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['42']
        self.A[ig,iv] = float(.826)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(14.61)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['70']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        ig = ig_['C']
        self.A[ig,iv] = float(.155)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['21']
        self.A[ig,iv] = float(-.3321)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['22']
        self.A[ig,iv] = float(-.5875)+self.A[ig,iv]
        ig = ig_['23']
        self.A[ig,iv] = float(-.362)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['39']
        self.A[ig,iv] = float(2.3)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.2049)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['42']
        self.A[ig,iv] = float(.826)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(14.61)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['66']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['70']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        ig = ig_['C']
        self.A[ig,iv] = float(.155)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['21']
        self.A[ig,iv] = float(-.2414)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['22']
        self.A[ig,iv] = float(-.6627)+self.A[ig,iv]
        ig = ig_['23']
        self.A[ig,iv] = float(-.293)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['39']
        self.A[ig,iv] = float(2.3)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.1531)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['42']
        self.A[ig,iv] = float(.826)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(14.61)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['70']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        ig = ig_['C']
        self.A[ig,iv] = float(.155)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        ig = ig_['21']
        self.A[ig,iv] = float(-.2414)+self.A[ig,iv]
        ig = ig_['22']
        self.A[ig,iv] = float(-.6627)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        ig = ig_['23']
        self.A[ig,iv] = float(-.293)+self.A[ig,iv]
        ig = ig_['28']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        ig = ig_['39']
        self.A[ig,iv] = float(2.3)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.1531)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        ig = ig_['42']
        self.A[ig,iv] = float(.826)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(14.61)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        ig = ig_['70']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.155)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.0185)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['13']
        self.A[ig,iv] = float(-.0568)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(-.0806)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['25']
        self.A[ig,iv] = float(-.0658)+self.A[ig,iv]
        ig = ig_['26']
        self.A[ig,iv] = float(-.0328)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['27']
        self.A[ig,iv] = float(-.4934)+self.A[ig,iv]
        ig = ig_['28']
        self.A[ig,iv] = float(-.2922)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['29']
        self.A[ig,iv] = float(-.0096)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(-.0654)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['41']
        self.A[ig,iv] = float(-.2535)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.632)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['43']
        self.A[ig,iv] = float(.6807)+self.A[ig,iv]
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        ig = ig_['71']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0528)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.0185)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['13']
        self.A[ig,iv] = float(-.0568)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(-.0806)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['25']
        self.A[ig,iv] = float(-.0658)+self.A[ig,iv]
        ig = ig_['26']
        self.A[ig,iv] = float(-.0328)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['27']
        self.A[ig,iv] = float(-.4934)+self.A[ig,iv]
        ig = ig_['28']
        self.A[ig,iv] = float(-.2922)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['29']
        self.A[ig,iv] = float(-.0096)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(-.0654)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['41']
        self.A[ig,iv] = float(-.2535)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(.632)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['43']
        self.A[ig,iv] = float(.6807)+self.A[ig,iv]
        ig = ig_['66']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        ig = ig_['71']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0528)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.0184)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['13']
        self.A[ig,iv] = float(-.0564)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(-.078)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['25']
        self.A[ig,iv] = float(-.0655)+self.A[ig,iv]
        ig = ig_['26']
        self.A[ig,iv] = float(-.0303)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['27']
        self.A[ig,iv] = float(-.475)+self.A[ig,iv]
        ig = ig_['28']
        self.A[ig,iv] = float(-.305)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['40']
        self.A[ig,iv] = float(-.0654)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.2703)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['42']
        self.A[ig,iv] = float(.632)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.6807)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['71']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        ig = ig_['C']
        self.A[ig,iv] = float(.0528)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.0184)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['13']
        self.A[ig,iv] = float(-.0564)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(-.078)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['25']
        self.A[ig,iv] = float(-.0655)+self.A[ig,iv]
        ig = ig_['26']
        self.A[ig,iv] = float(-.0303)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['27']
        self.A[ig,iv] = float(-.475)+self.A[ig,iv]
        ig = ig_['28']
        self.A[ig,iv] = float(-.305)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['40']
        self.A[ig,iv] = float(-.0654)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.2703)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['42']
        self.A[ig,iv] = float(.632)+self.A[ig,iv]
        ig = ig_['43']
        self.A[ig,iv] = float(.6807)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['71']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        ig = ig_['C']
        self.A[ig,iv] = float(.0528)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        ig = ig_['13']
        self.A[ig,iv] = float(.76)+self.A[ig,iv]
        ig = ig_['25']
        self.A[ig,iv] = float(.5714)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        ig = ig_['30']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.1869)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        ig = ig_['41']
        self.A[ig,iv] = float(.2796)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(2.241)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        ig = ig_['43']
        self.A[ig,iv] = float(2.766)+self.A[ig,iv]
        ig = ig_['72']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        ig = ig_['C']
        self.A[ig,iv] = float(.128)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['9']
        self.A[ig,iv] = float(-.0571)+self.A[ig,iv]
        ig = ig_['10']
        self.A[ig,iv] = float(-.0114)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['13']
        self.A[ig,iv] = float(.6571)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(.5714)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['31']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['40']
        self.A[ig,iv] = float(.1724)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['41']
        self.A[ig,iv] = float(.2579)+self.A[ig,iv]
        ig = ig_['42']
        self.A[ig,iv] = float(2.067)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['43']
        self.A[ig,iv] = float(2.552)+self.A[ig,iv]
        ig = ig_['72']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        ig = ig_['C']
        self.A[ig,iv] = float(.118)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('25',ix_)
        self.xnames=arrset(self.xnames,iv,'25')
        ig = ig_['9']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['25']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('26',ix_)
        self.xnames=arrset(self.xnames,iv,'26')
        ig = ig_['10']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['24']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('27',ix_)
        self.xnames=arrset(self.xnames,iv,'27')
        ig = ig_['10']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('28',ix_)
        self.xnames=arrset(self.xnames,iv,'28')
        ig = ig_['11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('28',ix_)
        self.xnames=arrset(self.xnames,iv,'28')
        ig = ig_['44']
        self.A[ig,iv] = float(-7.95)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-8.7)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('28',ix_)
        self.xnames=arrset(self.xnames,iv,'28')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(14.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('28',ix_)
        self.xnames=arrset(self.xnames,iv,'28')
        ig = ig_['48']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('29',ix_)
        self.xnames=arrset(self.xnames,iv,'29')
        ig = ig_['23']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('29',ix_)
        self.xnames=arrset(self.xnames,iv,'29')
        ig = ig_['44']
        self.A[ig,iv] = float(-8.84)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.45)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('29',ix_)
        self.xnames=arrset(self.xnames,iv,'29')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(12.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('29',ix_)
        self.xnames=arrset(self.xnames,iv,'29')
        ig = ig_['48']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('30',ix_)
        self.xnames=arrset(self.xnames,iv,'30')
        ig = ig_['19']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('30',ix_)
        self.xnames=arrset(self.xnames,iv,'30')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.43)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.57)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('30',ix_)
        self.xnames=arrset(self.xnames,iv,'30')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('30',ix_)
        self.xnames=arrset(self.xnames,iv,'30')
        ig = ig_['48']
        self.A[ig,iv] = float(.233)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-.358)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('31',ix_)
        self.xnames=arrset(self.xnames,iv,'31')
        ig = ig_['20']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('31',ix_)
        self.xnames=arrset(self.xnames,iv,'31')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.03)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.32)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('31',ix_)
        self.xnames=arrset(self.xnames,iv,'31')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('31',ix_)
        self.xnames=arrset(self.xnames,iv,'31')
        ig = ig_['48']
        self.A[ig,iv] = float(.205)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-.333)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('32',ix_)
        self.xnames=arrset(self.xnames,iv,'32')
        ig = ig_['27']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('32',ix_)
        self.xnames=arrset(self.xnames,iv,'32')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.23)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.22)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('32',ix_)
        self.xnames=arrset(self.xnames,iv,'32')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(6.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('32',ix_)
        self.xnames=arrset(self.xnames,iv,'32')
        ig = ig_['48']
        self.A[ig,iv] = float(.381)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-.509)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('33',ix_)
        self.xnames=arrset(self.xnames,iv,'33')
        ig = ig_['30']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('33',ix_)
        self.xnames=arrset(self.xnames,iv,'33')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.4)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.85)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('33',ix_)
        self.xnames=arrset(self.xnames,iv,'33')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(2.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('33',ix_)
        self.xnames=arrset(self.xnames,iv,'33')
        ig = ig_['48']
        self.A[ig,iv] = float(.39)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-.77)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('34',ix_)
        self.xnames=arrset(self.xnames,iv,'34')
        ig = ig_['31']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('34',ix_)
        self.xnames=arrset(self.xnames,iv,'34')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.74)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-10.1)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('34',ix_)
        self.xnames=arrset(self.xnames,iv,'34')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(3.3)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('34',ix_)
        self.xnames=arrset(self.xnames,iv,'34')
        ig = ig_['48']
        self.A[ig,iv] = float(.233)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-.58)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('35',ix_)
        self.xnames=arrset(self.xnames,iv,'35')
        ig = ig_['10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('35',ix_)
        self.xnames=arrset(self.xnames,iv,'35')
        ig = ig_['44']
        self.A[ig,iv] = float(-9.74)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-9.9)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('35',ix_)
        self.xnames=arrset(self.xnames,iv,'35')
        ig = ig_['46']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(66.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('35',ix_)
        self.xnames=arrset(self.xnames,iv,'35')
        ig = ig_['48']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('36',ix_)
        self.xnames=arrset(self.xnames,iv,'36')
        ig = ig_['44']
        self.A[ig,iv] = float(-.493)+self.A[ig,iv]
        ig = ig_['45']
        self.A[ig,iv] = float(-.165)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('36',ix_)
        self.xnames=arrset(self.xnames,iv,'36')
        ig = ig_['46']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0924)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        ig = ig_['32']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['44']
        self.A[ig,iv] = float(10.03)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        ig = ig_['45']
        self.A[ig,iv] = float(10.03)+self.A[ig,iv]
        ig = ig_['47']
        self.A[ig,iv] = float(-9.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        ig = ig_['48']
        self.A[ig,iv] = float(-.5)+self.A[ig,iv]
        ig = ig_['49']
        self.A[ig,iv] = float(.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        ig = ig_['73']
        self.A[ig,iv] = float(.64)+self.A[ig,iv]
        ig = ig_['74']
        self.A[ig,iv] = float(.35)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        ig = ig_['C']
        self.A[ig,iv] = float(-5.36)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('38',ix_)
        self.xnames=arrset(self.xnames,iv,'38')
        ig = ig_['11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('38',ix_)
        self.xnames=arrset(self.xnames,iv,'38')
        ig = ig_['50']
        self.A[ig,iv] = float(-7.98)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-8.58)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('38',ix_)
        self.xnames=arrset(self.xnames,iv,'38')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(14.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('38',ix_)
        self.xnames=arrset(self.xnames,iv,'38')
        ig = ig_['54']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('39',ix_)
        self.xnames=arrset(self.xnames,iv,'39')
        ig = ig_['23']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('39',ix_)
        self.xnames=arrset(self.xnames,iv,'39')
        ig = ig_['50']
        self.A[ig,iv] = float(-8.87)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-9.33)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('39',ix_)
        self.xnames=arrset(self.xnames,iv,'39')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(12.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('39',ix_)
        self.xnames=arrset(self.xnames,iv,'39')
        ig = ig_['54']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('40',ix_)
        self.xnames=arrset(self.xnames,iv,'40')
        ig = ig_['19']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('40',ix_)
        self.xnames=arrset(self.xnames,iv,'40')
        ig = ig_['50']
        self.A[ig,iv] = float(-9.46)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-9.45)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('40',ix_)
        self.xnames=arrset(self.xnames,iv,'40')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('40',ix_)
        self.xnames=arrset(self.xnames,iv,'40')
        ig = ig_['54']
        self.A[ig,iv] = float(.233)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-.358)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('41',ix_)
        self.xnames=arrset(self.xnames,iv,'41')
        ig = ig_['20']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('41',ix_)
        self.xnames=arrset(self.xnames,iv,'41')
        ig = ig_['50']
        self.A[ig,iv] = float(-9.06)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-9.2)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('41',ix_)
        self.xnames=arrset(self.xnames,iv,'41')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('41',ix_)
        self.xnames=arrset(self.xnames,iv,'41')
        ig = ig_['54']
        self.A[ig,iv] = float(.205)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-.333)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('42',ix_)
        self.xnames=arrset(self.xnames,iv,'42')
        ig = ig_['27']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('42',ix_)
        self.xnames=arrset(self.xnames,iv,'42')
        ig = ig_['50']
        self.A[ig,iv] = float(-9.26)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-9.13)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('42',ix_)
        self.xnames=arrset(self.xnames,iv,'42')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(6.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('42',ix_)
        self.xnames=arrset(self.xnames,iv,'42')
        ig = ig_['54']
        self.A[ig,iv] = float(.318)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-.509)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('43',ix_)
        self.xnames=arrset(self.xnames,iv,'43')
        ig = ig_['10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('43',ix_)
        self.xnames=arrset(self.xnames,iv,'43')
        ig = ig_['50']
        self.A[ig,iv] = float(-9.77)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-9.78)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('43',ix_)
        self.xnames=arrset(self.xnames,iv,'43')
        ig = ig_['52']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(66.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('43',ix_)
        self.xnames=arrset(self.xnames,iv,'43')
        ig = ig_['54']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('44',ix_)
        self.xnames=arrset(self.xnames,iv,'44')
        ig = ig_['50']
        self.A[ig,iv] = float(-.435)+self.A[ig,iv]
        ig = ig_['51']
        self.A[ig,iv] = float(-.208)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('44',ix_)
        self.xnames=arrset(self.xnames,iv,'44')
        ig = ig_['52']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0924)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        ig = ig_['33']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['50']
        self.A[ig,iv] = float(9.65)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        ig = ig_['51']
        self.A[ig,iv] = float(9.65)+self.A[ig,iv]
        ig = ig_['53']
        self.A[ig,iv] = float(-9.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        ig = ig_['54']
        self.A[ig,iv] = float(-.5)+self.A[ig,iv]
        ig = ig_['55']
        self.A[ig,iv] = float(.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        ig = ig_['73']
        self.A[ig,iv] = float(-.36)+self.A[ig,iv]
        ig = ig_['74']
        self.A[ig,iv] = float(.35)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        ig = ig_['C']
        self.A[ig,iv] = float(-5.08)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('46',ix_)
        self.xnames=arrset(self.xnames,iv,'46')
        ig = ig_['11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('46',ix_)
        self.xnames=arrset(self.xnames,iv,'46')
        ig = ig_['56']
        self.A[ig,iv] = float(-7.99)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-8.59)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('46',ix_)
        self.xnames=arrset(self.xnames,iv,'46')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(14.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('46',ix_)
        self.xnames=arrset(self.xnames,iv,'46')
        ig = ig_['60']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('47',ix_)
        self.xnames=arrset(self.xnames,iv,'47')
        ig = ig_['23']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('47',ix_)
        self.xnames=arrset(self.xnames,iv,'47')
        ig = ig_['56']
        self.A[ig,iv] = float(-8.88)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-9.34)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('47',ix_)
        self.xnames=arrset(self.xnames,iv,'47')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(12.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('47',ix_)
        self.xnames=arrset(self.xnames,iv,'47')
        ig = ig_['60']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('48',ix_)
        self.xnames=arrset(self.xnames,iv,'48')
        ig = ig_['19']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('48',ix_)
        self.xnames=arrset(self.xnames,iv,'48')
        ig = ig_['56']
        self.A[ig,iv] = float(-9.47)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-9.46)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('48',ix_)
        self.xnames=arrset(self.xnames,iv,'48')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('48',ix_)
        self.xnames=arrset(self.xnames,iv,'48')
        ig = ig_['60']
        self.A[ig,iv] = float(.233)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-.358)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('49',ix_)
        self.xnames=arrset(self.xnames,iv,'49')
        ig = ig_['20']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('49',ix_)
        self.xnames=arrset(self.xnames,iv,'49')
        ig = ig_['56']
        self.A[ig,iv] = float(-9.07)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-9.21)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('49',ix_)
        self.xnames=arrset(self.xnames,iv,'49')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(3.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('49',ix_)
        self.xnames=arrset(self.xnames,iv,'49')
        ig = ig_['60']
        self.A[ig,iv] = float(.205)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-.333)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('50',ix_)
        self.xnames=arrset(self.xnames,iv,'50')
        ig = ig_['27']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('50',ix_)
        self.xnames=arrset(self.xnames,iv,'50')
        ig = ig_['56']
        self.A[ig,iv] = float(-9.27)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-9.14)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('50',ix_)
        self.xnames=arrset(self.xnames,iv,'50')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(6.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('50',ix_)
        self.xnames=arrset(self.xnames,iv,'50')
        ig = ig_['60']
        self.A[ig,iv] = float(.318)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-.509)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('51',ix_)
        self.xnames=arrset(self.xnames,iv,'51')
        ig = ig_['10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('51',ix_)
        self.xnames=arrset(self.xnames,iv,'51')
        ig = ig_['56']
        self.A[ig,iv] = float(-9.78)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-9.79)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('51',ix_)
        self.xnames=arrset(self.xnames,iv,'51')
        ig = ig_['58']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(66.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('51',ix_)
        self.xnames=arrset(self.xnames,iv,'51')
        ig = ig_['60']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('52',ix_)
        self.xnames=arrset(self.xnames,iv,'52')
        ig = ig_['56']
        self.A[ig,iv] = float(-.426)+self.A[ig,iv]
        ig = ig_['57']
        self.A[ig,iv] = float(-.204)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('52',ix_)
        self.xnames=arrset(self.xnames,iv,'52')
        ig = ig_['58']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0924)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('53',ix_)
        self.xnames=arrset(self.xnames,iv,'53')
        ig = ig_['36']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['56']
        self.A[ig,iv] = float(9.05)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('53',ix_)
        self.xnames=arrset(self.xnames,iv,'53')
        ig = ig_['57']
        self.A[ig,iv] = float(9.05)+self.A[ig,iv]
        ig = ig_['59']
        self.A[ig,iv] = float(-9.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('53',ix_)
        self.xnames=arrset(self.xnames,iv,'53')
        ig = ig_['60']
        self.A[ig,iv] = float(-.5)+self.A[ig,iv]
        ig = ig_['61']
        self.A[ig,iv] = float(.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('53',ix_)
        self.xnames=arrset(self.xnames,iv,'53')
        ig = ig_['73']
        self.A[ig,iv] = float(-.36)+self.A[ig,iv]
        ig = ig_['74']
        self.A[ig,iv] = float(-.65)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('53',ix_)
        self.xnames=arrset(self.xnames,iv,'53')
        ig = ig_['C']
        self.A[ig,iv] = float(-4.51)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('54',ix_)
        self.xnames=arrset(self.xnames,iv,'54')
        ig = ig_['9']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['26']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('55',ix_)
        self.xnames=arrset(self.xnames,iv,'55')
        ig = ig_['9']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['37']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('56',ix_)
        self.xnames=arrset(self.xnames,iv,'56')
        ig = ig_['10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['37']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('57',ix_)
        self.xnames=arrset(self.xnames,iv,'57')
        ig = ig_['37']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(-2.75)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('58',ix_)
        self.xnames=arrset(self.xnames,iv,'58')
        ig = ig_['11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['38']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('58',ix_)
        self.xnames=arrset(self.xnames,iv,'58')
        ig = ig_['63']
        self.A[ig,iv] = float(-14.0)+self.A[ig,iv]
        ig = ig_['64']
        self.A[ig,iv] = float(14.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('59',ix_)
        self.xnames=arrset(self.xnames,iv,'59')
        ig = ig_['12']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['38']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('59',ix_)
        self.xnames=arrset(self.xnames,iv,'59')
        ig = ig_['63']
        self.A[ig,iv] = float(-.8)+self.A[ig,iv]
        ig = ig_['64']
        self.A[ig,iv] = float(.8)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('60',ix_)
        self.xnames=arrset(self.xnames,iv,'60')
        ig = ig_['38']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['63']
        self.A[ig,iv] = float(2.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('60',ix_)
        self.xnames=arrset(self.xnames,iv,'60')
        ig = ig_['64']
        self.A[ig,iv] = float(-3.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(-4.2)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('61',ix_)
        self.xnames=arrset(self.xnames,iv,'61')
        ig = ig_['4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['34']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('62',ix_)
        self.xnames=arrset(self.xnames,iv,'62')
        ig = ig_['3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['34']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('63',ix_)
        self.xnames=arrset(self.xnames,iv,'63')
        ig = ig_['34']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('63',ix_)
        self.xnames=arrset(self.xnames,iv,'63')
        ig = ig_['C']
        self.A[ig,iv] = float(-3.6)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('64',ix_)
        self.xnames=arrset(self.xnames,iv,'64')
        ig = ig_['7']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('64',ix_)
        self.xnames=arrset(self.xnames,iv,'64')
        ig = ig_['62']
        self.A[ig,iv] = float(10.1)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('65',ix_)
        self.xnames=arrset(self.xnames,iv,'65')
        ig = ig_['8']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('65',ix_)
        self.xnames=arrset(self.xnames,iv,'65')
        ig = ig_['62']
        self.A[ig,iv] = float(12.63)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('66',ix_)
        self.xnames=arrset(self.xnames,iv,'66')
        ig = ig_['6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('66',ix_)
        self.xnames=arrset(self.xnames,iv,'66')
        ig = ig_['62']
        self.A[ig,iv] = float(8.05)+self.A[ig,iv]
        ig = ig_['66']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('67',ix_)
        self.xnames=arrset(self.xnames,iv,'67')
        ig = ig_['5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('67',ix_)
        self.xnames=arrset(self.xnames,iv,'67')
        ig = ig_['62']
        self.A[ig,iv] = float(6.9)+self.A[ig,iv]
        ig = ig_['65']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('68',ix_)
        self.xnames=arrset(self.xnames,iv,'68')
        ig = ig_['29']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('68',ix_)
        self.xnames=arrset(self.xnames,iv,'68')
        ig = ig_['62']
        self.A[ig,iv] = float(8.05)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('69',ix_)
        self.xnames=arrset(self.xnames,iv,'69')
        ig = ig_['28']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('69',ix_)
        self.xnames=arrset(self.xnames,iv,'69')
        ig = ig_['62']
        self.A[ig,iv] = float(4.4)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('70',ix_)
        self.xnames=arrset(self.xnames,iv,'70')
        ig = ig_['35']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['62']
        self.A[ig,iv] = float(-10.1)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('70',ix_)
        self.xnames=arrset(self.xnames,iv,'70')
        ig = ig_['C']
        self.A[ig,iv] = float(-2.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('71',ix_)
        self.xnames=arrset(self.xnames,iv,'71')
        ig = ig_['39']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-.325)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('72',ix_)
        self.xnames=arrset(self.xnames,iv,'72')
        ig = ig_['13']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-4.153)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('73',ix_)
        self.xnames=arrset(self.xnames,iv,'73')
        ig = ig_['10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-4.316)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('74',ix_)
        self.xnames=arrset(self.xnames,iv,'74')
        ig = ig_['9']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-3.814)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('75',ix_)
        self.xnames=arrset(self.xnames,iv,'75')
        ig = ig_['25']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-3.808)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('76',ix_)
        self.xnames=arrset(self.xnames,iv,'76')
        ig = ig_['24']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(-4.44)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('77',ix_)
        self.xnames=arrset(self.xnames,iv,'77')
        ig = ig_['40']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['41']
        self.A[ig,iv] = float(1.42)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('77',ix_)
        self.xnames=arrset(self.xnames,iv,'77')
        ig = ig_['C']
        self.A[ig,iv] = float(.04)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('78',ix_)
        self.xnames=arrset(self.xnames,iv,'78')
        ig = ig_['40']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('79',ix_)
        self.xnames=arrset(self.xnames,iv,'79')
        ig = ig_['10']
        self.A[ig,iv] = float(-.5)+self.A[ig,iv]
        ig = ig_['13']
        self.A[ig,iv] = float(-.5)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('79',ix_)
        self.xnames=arrset(self.xnames,iv,'79')
        ig = ig_['C']
        self.A[ig,iv] = float(3.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('80',ix_)
        self.xnames=arrset(self.xnames,iv,'80')
        ig = ig_['41']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.4)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('81',ix_)
        self.xnames=arrset(self.xnames,iv,'81')
        ig = ig_['41']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('82',ix_)
        self.xnames=arrset(self.xnames,iv,'82')
        ig = ig_['42']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.0132)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('83',ix_)
        self.xnames=arrset(self.xnames,iv,'83')
        ig = ig_['43']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        ig = ig_['C']
        self.A[ig,iv] = float(.01)+self.A[ig,iv]
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
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['65'],float(23.26))
        self.gconst = arrset(self.gconst,ig_['66'],float(5.25))
        self.gconst = arrset(self.gconst,ig_['67'],float(26.32))
        self.gconst = arrset(self.gconst,ig_['68'],float(21.05))
        self.gconst = arrset(self.gconst,ig_['69'],float(13.45))
        self.gconst = arrset(self.gconst,ig_['70'],float(2.58))
        self.gconst = arrset(self.gconst,ig_['71'],float(10.0))
        self.gconst = arrset(self.gconst,ig_['72'],float(10.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'D')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        self.x0 = np.zeros((self.n,1))
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.000000000))
        ename = 'E2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.109756112))
        ename = 'E3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.219512224))
        ename = 'E4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.329268336))
        ename = 'E5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.439024448))
        ename = 'E6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.548780441))
        ename = 'E7'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.658536553))
        ename = 'E8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.768292665))
        ename = 'E9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.878048778))
        ename = 'E10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(1.987804890))
        ename = 'E11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.097560883))
        ename = 'E12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.207317114))
        ename = 'E13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.317073107))
        ename = 'E14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.426829338))
        ename = 'E15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.536585331))
        ename = 'E16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.646341562))
        ename = 'E17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.756097555))
        ename = 'E18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '18'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.865853548))
        ename = 'E19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '19'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(2.975609779))
        ename = 'E20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.085365772))
        ename = 'E21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '21'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.195122004))
        ename = 'E22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '22'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.304877996))
        ename = 'E23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '23'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.414634228))
        ename = 'E24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '24'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.524390221))
        ename = 'E25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '25'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.634146452))
        ename = 'E26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '26'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.743902445))
        ename = 'E27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '27'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.853658438))
        ename = 'E28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '28'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(3.963414669))
        ename = 'E29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.073170662))
        ename = 'E30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '30'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.182926655))
        ename = 'E31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '31'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.292683125))
        ename = 'E32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.402439117))
        ename = 'E33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '33'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.512195110))
        ename = 'E34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '34'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.621951103))
        ename = 'E35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '35'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.731707096))
        ename = 'E36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.841463566))
        ename = 'E37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '37'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(4.951219559))
        ename = 'E38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '38'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.060975552))
        ename = 'E39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '39'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.170731544))
        ename = 'E40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '40'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.280488014))
        ename = 'E41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.390244007))
        ename = 'E42'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '42'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.500000000))
        ename = 'E43'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '43'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.609755993))
        ename = 'E44'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '44'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.719511986))
        ename = 'E45'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '45'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.829268456))
        ename = 'E46'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '46'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(5.939024448))
        ename = 'E47'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '47'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.048780441))
        ename = 'E48'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '48'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.158536434))
        ename = 'E49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.268292904))
        ename = 'E50'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '50'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.378048897))
        ename = 'E51'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '51'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.487804890))
        ename = 'E52'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '52'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.597560883))
        ename = 'E53'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '53'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.707316875))
        ename = 'E54'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '54'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.817073345))
        ename = 'E55'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '55'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(6.926829338))
        ename = 'E56'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '56'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.036585331))
        ename = 'E57'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '57'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.146341324))
        ename = 'E58'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '58'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.256097794))
        ename = 'E59'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '59'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.365853786))
        ename = 'E60'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '60'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.475609779))
        ename = 'E61'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '61'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.585365772))
        ename = 'E62'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '62'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.695121765))
        ename = 'E63'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '63'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.804878235))
        ename = 'E64'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '64'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(7.914634228))
        ename = 'E65'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '65'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.024390221))
        ename = 'E66'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '66'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.134146690))
        ename = 'E67'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '67'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.243902206))
        ename = 'E68'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '68'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.353658676))
        ename = 'E69'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '69'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.463414192))
        ename = 'E70'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '70'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.573170662))
        ename = 'E71'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '71'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.682927132))
        ename = 'E72'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '72'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.792682648))
        ename = 'E73'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '73'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(8.902439117))
        ename = 'E74'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '74'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.012195587))
        ename = 'E75'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '75'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.121951103))
        ename = 'E76'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '76'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.231707573))
        ename = 'E77'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '77'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.341463089))
        ename = 'E78'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '78'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.451219559))
        ename = 'E79'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '79'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.560976028))
        ename = 'E80'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '80'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.670731544))
        ename = 'E81'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '81'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.780488014))
        ename = 'E82'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '82'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(9.890243530))
        ename = 'E83'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eSQR')
            ielftype = arrset(ielftype,ie,iet_['eSQR'])
        vname = '83'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='D')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(10.000000000))
        v_['1'] = 1
        v_['N'] = 83
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['C']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass = "C-CQLR2-MN-83-74"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*self.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0*self.elpar[iel_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

