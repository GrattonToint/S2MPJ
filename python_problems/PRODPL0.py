from s2mpjlib import *
class  PRODPL0(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PRODPL0
#    *********
# 
#    A production planning problem in the computer industry.
# 
#    Source:
#    L. Escudero, private communication, 1991.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-CLQR2-RY-60-29"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PRODPL0'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['T'] = 5
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('COST',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('K01',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K01')
        [ig,ig_,_] = s2mpj_ii('K02',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K02')
        [ig,ig_,_] = s2mpj_ii('K03',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K03')
        [ig,ig_,_] = s2mpj_ii('K04',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K04')
        [ig,ig_,_] = s2mpj_ii('K05',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K05')
        [ig,ig_,_] = s2mpj_ii('D00101',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00101')
        [ig,ig_,_] = s2mpj_ii('D00201',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00201')
        [ig,ig_,_] = s2mpj_ii('D00301',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00301')
        [ig,ig_,_] = s2mpj_ii('D00401',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00401')
        [ig,ig_,_] = s2mpj_ii('D00102',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00102')
        [ig,ig_,_] = s2mpj_ii('D00202',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00202')
        [ig,ig_,_] = s2mpj_ii('D00302',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00302')
        [ig,ig_,_] = s2mpj_ii('D00402',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00402')
        [ig,ig_,_] = s2mpj_ii('D00103',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00103')
        [ig,ig_,_] = s2mpj_ii('D00203',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00203')
        [ig,ig_,_] = s2mpj_ii('D00303',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00303')
        [ig,ig_,_] = s2mpj_ii('D00403',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00403')
        [ig,ig_,_] = s2mpj_ii('D00104',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00104')
        [ig,ig_,_] = s2mpj_ii('D00204',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00204')
        [ig,ig_,_] = s2mpj_ii('D00304',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00304')
        [ig,ig_,_] = s2mpj_ii('D00404',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00404')
        [ig,ig_,_] = s2mpj_ii('D00105',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00105')
        [ig,ig_,_] = s2mpj_ii('D00205',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00205')
        [ig,ig_,_] = s2mpj_ii('D00305',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00305')
        [ig,ig_,_] = s2mpj_ii('D00405',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00405')
        v_['TM1'] = -1+v_['T']
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            [ig,ig_,_] = s2mpj_ii('SMOOTH'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SMOOTH'+str(I))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('X00101',ix_)
        self.xnames=arrset(self.xnames,iv,'X00101')
        ig = ig_['K01']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00101']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00201',ix_)
        self.xnames=arrset(self.xnames,iv,'X00201')
        ig = ig_['K01']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00201']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00301',ix_)
        self.xnames=arrset(self.xnames,iv,'X00301')
        ig = ig_['K01']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00301']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00401',ix_)
        self.xnames=arrset(self.xnames,iv,'X00401')
        ig = ig_['K01']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00401']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00102',ix_)
        self.xnames=arrset(self.xnames,iv,'X00102')
        ig = ig_['K02']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00102']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00202',ix_)
        self.xnames=arrset(self.xnames,iv,'X00202')
        ig = ig_['K02']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00202']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00302',ix_)
        self.xnames=arrset(self.xnames,iv,'X00302')
        ig = ig_['K02']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00302']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00402',ix_)
        self.xnames=arrset(self.xnames,iv,'X00402')
        ig = ig_['K02']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00402']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00103',ix_)
        self.xnames=arrset(self.xnames,iv,'X00103')
        ig = ig_['K03']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00103']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00203',ix_)
        self.xnames=arrset(self.xnames,iv,'X00203')
        ig = ig_['K03']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00203']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00303',ix_)
        self.xnames=arrset(self.xnames,iv,'X00303')
        ig = ig_['K03']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00303']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00403',ix_)
        self.xnames=arrset(self.xnames,iv,'X00403')
        ig = ig_['K03']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00403']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00104',ix_)
        self.xnames=arrset(self.xnames,iv,'X00104')
        ig = ig_['K04']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00104']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00204',ix_)
        self.xnames=arrset(self.xnames,iv,'X00204')
        ig = ig_['K04']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00204']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00304',ix_)
        self.xnames=arrset(self.xnames,iv,'X00304')
        ig = ig_['K04']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00304']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00404',ix_)
        self.xnames=arrset(self.xnames,iv,'X00404')
        ig = ig_['K04']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00404']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00105',ix_)
        self.xnames=arrset(self.xnames,iv,'X00105')
        ig = ig_['K05']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00105']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00205',ix_)
        self.xnames=arrset(self.xnames,iv,'X00205')
        ig = ig_['K05']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00205']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00305',ix_)
        self.xnames=arrset(self.xnames,iv,'X00305')
        ig = ig_['K05']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00305']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('X00405',ix_)
        self.xnames=arrset(self.xnames,iv,'X00405')
        ig = ig_['K05']
        self.A[ig,iv] = float(1)+self.A[ig,iv]
        ig = ig_['D00405']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00101',ix_)
        self.xnames=arrset(self.xnames,iv,'I00101')
        ig = ig_['COST']
        self.A[ig,iv] = float(1.000000)+self.A[ig,iv]
        ig = ig_['D00101']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00101',ix_)
        self.xnames=arrset(self.xnames,iv,'I00101')
        ig = ig_['D00102']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00101',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00101')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00101']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00201',ix_)
        self.xnames=arrset(self.xnames,iv,'I00201')
        ig = ig_['COST']
        self.A[ig,iv] = float(2.000000)+self.A[ig,iv]
        ig = ig_['D00201']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00201',ix_)
        self.xnames=arrset(self.xnames,iv,'I00201')
        ig = ig_['D00202']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00201',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00201')
        ig = ig_['COST']
        self.A[ig,iv] = float(3)+self.A[ig,iv]
        ig = ig_['D00201']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00301',ix_)
        self.xnames=arrset(self.xnames,iv,'I00301')
        ig = ig_['COST']
        self.A[ig,iv] = float(3.000000)+self.A[ig,iv]
        ig = ig_['D00301']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00301',ix_)
        self.xnames=arrset(self.xnames,iv,'I00301')
        ig = ig_['D00302']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00301',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00301')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00301']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00401',ix_)
        self.xnames=arrset(self.xnames,iv,'I00401')
        ig = ig_['COST']
        self.A[ig,iv] = float(4.000000)+self.A[ig,iv]
        ig = ig_['D00401']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00401',ix_)
        self.xnames=arrset(self.xnames,iv,'I00401')
        ig = ig_['D00402']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00401',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00401')
        ig = ig_['COST']
        self.A[ig,iv] = float(5)+self.A[ig,iv]
        ig = ig_['D00401']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00102',ix_)
        self.xnames=arrset(self.xnames,iv,'I00102')
        ig = ig_['COST']
        self.A[ig,iv] = float(1.000000)+self.A[ig,iv]
        ig = ig_['D00102']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00102',ix_)
        self.xnames=arrset(self.xnames,iv,'I00102')
        ig = ig_['D00103']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00102',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00102')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00102']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00202',ix_)
        self.xnames=arrset(self.xnames,iv,'I00202')
        ig = ig_['COST']
        self.A[ig,iv] = float(2.000000)+self.A[ig,iv]
        ig = ig_['D00202']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00202',ix_)
        self.xnames=arrset(self.xnames,iv,'I00202')
        ig = ig_['D00203']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00202',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00202')
        ig = ig_['COST']
        self.A[ig,iv] = float(3)+self.A[ig,iv]
        ig = ig_['D00202']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00302',ix_)
        self.xnames=arrset(self.xnames,iv,'I00302')
        ig = ig_['COST']
        self.A[ig,iv] = float(3.000000)+self.A[ig,iv]
        ig = ig_['D00302']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00302',ix_)
        self.xnames=arrset(self.xnames,iv,'I00302')
        ig = ig_['D00303']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00302',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00302')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00302']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00402',ix_)
        self.xnames=arrset(self.xnames,iv,'I00402')
        ig = ig_['COST']
        self.A[ig,iv] = float(4.000000)+self.A[ig,iv]
        ig = ig_['D00402']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00402',ix_)
        self.xnames=arrset(self.xnames,iv,'I00402')
        ig = ig_['D00403']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00402',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00402')
        ig = ig_['COST']
        self.A[ig,iv] = float(5)+self.A[ig,iv]
        ig = ig_['D00402']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00103',ix_)
        self.xnames=arrset(self.xnames,iv,'I00103')
        ig = ig_['COST']
        self.A[ig,iv] = float(1.000000)+self.A[ig,iv]
        ig = ig_['D00103']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00103',ix_)
        self.xnames=arrset(self.xnames,iv,'I00103')
        ig = ig_['D00104']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00103',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00103')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00103']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00203',ix_)
        self.xnames=arrset(self.xnames,iv,'I00203')
        ig = ig_['COST']
        self.A[ig,iv] = float(2.000000)+self.A[ig,iv]
        ig = ig_['D00203']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00203',ix_)
        self.xnames=arrset(self.xnames,iv,'I00203')
        ig = ig_['D00204']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00203',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00203')
        ig = ig_['COST']
        self.A[ig,iv] = float(3)+self.A[ig,iv]
        ig = ig_['D00203']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00303',ix_)
        self.xnames=arrset(self.xnames,iv,'I00303')
        ig = ig_['COST']
        self.A[ig,iv] = float(3.000000)+self.A[ig,iv]
        ig = ig_['D00303']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00303',ix_)
        self.xnames=arrset(self.xnames,iv,'I00303')
        ig = ig_['D00304']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00303',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00303')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00303']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00403',ix_)
        self.xnames=arrset(self.xnames,iv,'I00403')
        ig = ig_['COST']
        self.A[ig,iv] = float(4.000000)+self.A[ig,iv]
        ig = ig_['D00403']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00403',ix_)
        self.xnames=arrset(self.xnames,iv,'I00403')
        ig = ig_['D00404']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00403',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00403')
        ig = ig_['COST']
        self.A[ig,iv] = float(5)+self.A[ig,iv]
        ig = ig_['D00403']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00104',ix_)
        self.xnames=arrset(self.xnames,iv,'I00104')
        ig = ig_['COST']
        self.A[ig,iv] = float(1.000000)+self.A[ig,iv]
        ig = ig_['D00104']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00104',ix_)
        self.xnames=arrset(self.xnames,iv,'I00104')
        ig = ig_['D00105']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00104',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00104')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00104']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00204',ix_)
        self.xnames=arrset(self.xnames,iv,'I00204')
        ig = ig_['COST']
        self.A[ig,iv] = float(2.000000)+self.A[ig,iv]
        ig = ig_['D00204']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00204',ix_)
        self.xnames=arrset(self.xnames,iv,'I00204')
        ig = ig_['D00205']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00204',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00204')
        ig = ig_['COST']
        self.A[ig,iv] = float(3)+self.A[ig,iv]
        ig = ig_['D00204']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00304',ix_)
        self.xnames=arrset(self.xnames,iv,'I00304')
        ig = ig_['COST']
        self.A[ig,iv] = float(3.000000)+self.A[ig,iv]
        ig = ig_['D00304']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00304',ix_)
        self.xnames=arrset(self.xnames,iv,'I00304')
        ig = ig_['D00305']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00304',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00304')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00304']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00404',ix_)
        self.xnames=arrset(self.xnames,iv,'I00404')
        ig = ig_['COST']
        self.A[ig,iv] = float(4.000000)+self.A[ig,iv]
        ig = ig_['D00404']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00404',ix_)
        self.xnames=arrset(self.xnames,iv,'I00404')
        ig = ig_['D00405']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00404',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00404')
        ig = ig_['COST']
        self.A[ig,iv] = float(5)+self.A[ig,iv]
        ig = ig_['D00404']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00105',ix_)
        self.xnames=arrset(self.xnames,iv,'I00105')
        ig = ig_['COST']
        self.A[ig,iv] = float(1.000000)+self.A[ig,iv]
        ig = ig_['D00105']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00105',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00105')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00105']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00205',ix_)
        self.xnames=arrset(self.xnames,iv,'I00205')
        ig = ig_['COST']
        self.A[ig,iv] = float(2.000000)+self.A[ig,iv]
        ig = ig_['D00205']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00205',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00205')
        ig = ig_['COST']
        self.A[ig,iv] = float(3)+self.A[ig,iv]
        ig = ig_['D00205']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00305',ix_)
        self.xnames=arrset(self.xnames,iv,'I00305')
        ig = ig_['COST']
        self.A[ig,iv] = float(3.000000)+self.A[ig,iv]
        ig = ig_['D00305']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00305',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00305')
        ig = ig_['COST']
        self.A[ig,iv] = float(2)+self.A[ig,iv]
        ig = ig_['D00305']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('I00405',ix_)
        self.xnames=arrset(self.xnames,iv,'I00405')
        ig = ig_['COST']
        self.A[ig,iv] = float(4.000000)+self.A[ig,iv]
        ig = ig_['D00405']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [iv,ix_,_] = s2mpj_ii('Y00405',ix_)
        self.xnames=arrset(self.xnames,iv,'Y00405')
        ig = ig_['COST']
        self.A[ig,iv] = float(5)+self.A[ig,iv]
        ig = ig_['D00405']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
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
        self.gconst = arrset(self.gconst,ig_['K01'],float(3))
        self.gconst = arrset(self.gconst,ig_['K02'],float(6))
        self.gconst = arrset(self.gconst,ig_['K03'],float(10))
        self.gconst = arrset(self.gconst,ig_['K04'],float(2000))
        self.gconst = arrset(self.gconst,ig_['K05'],float(18))
        self.gconst = arrset(self.gconst,ig_['D00101'],float(1.000))
        self.gconst = arrset(self.gconst,ig_['D00201'],float(1.000))
        self.gconst = arrset(self.gconst,ig_['D00301'],float(1.000))
        self.gconst = arrset(self.gconst,ig_['D00401'],float(1.000))
        self.gconst = arrset(self.gconst,ig_['D00102'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00202'],float(1.667))
        self.gconst = arrset(self.gconst,ig_['D00302'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00402'],float(3.333))
        self.gconst = arrset(self.gconst,ig_['D00103'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00203'],float(2.000))
        self.gconst = arrset(self.gconst,ig_['D00303'],float(3.000))
        self.gconst = arrset(self.gconst,ig_['D00403'],float(3.000))
        self.gconst = arrset(self.gconst,ig_['D00104'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00204'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00304'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00404'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00105'],float(2.667))
        self.gconst = arrset(self.gconst,ig_['D00205'],float(2.333))
        self.gconst = arrset(self.gconst,ig_['D00305'],float(2.333))
        self.gconst = arrset(self.gconst,ig_['D00405'],float(2.333))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQMRSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            v_['IP1'] = 1+I
            ename = 'NLE'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQMRSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQMRSQ"])
            self.x0 = np.zeros((self.n,1))
            vname = 'X0010'+str(int(v_['IP1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0020'+str(int(v_['IP1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0030'+str(int(v_['IP1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0040'+str(int(v_['IP1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0010'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0020'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0030'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V7')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X0040'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            ig = ig_['SMOOTH'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['NLE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               58.7898356794
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
        self.pbclass = "C-CLQR2-RY-60-29"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,2.0)
        self.efpar = arrset( self.efpar,1,0.1)
        self.efpar = arrset( self.efpar,2,self.efpar[1]*self.efpar[1])
        return pbm

    @staticmethod
    def eSQMRSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,8))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        U_[0,3] = U_[0,3]+1
        U_[1,4] = U_[1,4]+1
        U_[1,5] = U_[1,5]+1
        U_[1,6] = U_[1,6]+1
        U_[1,7] = U_[1,7]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        U1MU2 = IV_[0]-IV_[1]
        f_   = U1MU2**2-self.efpar[2]*IV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.efpar[0]*U1MU2
            g_[1] = -self.efpar[0]*U1MU2-self.efpar[0]*self.efpar[2]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = self.efpar[0]
                H_[0,1] = -self.efpar[0]
                H_[1,0] = H_[0,1]
                H_[1,1] = self.efpar[0]*(1.0-self.efpar[2])
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

