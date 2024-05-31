from s2xlib import *
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
#    classification = "LQR2-RY-60-29"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PRODPL0'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'PRODPL0'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['T'] = 5
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('COST',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2x_ii('K01',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K01')
        [ig,ig_,_] = s2x_ii('K02',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K02')
        [ig,ig_,_] = s2x_ii('K03',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K03')
        [ig,ig_,_] = s2x_ii('K04',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K04')
        [ig,ig_,_] = s2x_ii('K05',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'K05')
        [ig,ig_,_] = s2x_ii('D00101',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00101')
        [ig,ig_,_] = s2x_ii('D00201',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00201')
        [ig,ig_,_] = s2x_ii('D00301',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00301')
        [ig,ig_,_] = s2x_ii('D00401',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00401')
        [ig,ig_,_] = s2x_ii('D00102',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00102')
        [ig,ig_,_] = s2x_ii('D00202',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00202')
        [ig,ig_,_] = s2x_ii('D00302',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00302')
        [ig,ig_,_] = s2x_ii('D00402',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00402')
        [ig,ig_,_] = s2x_ii('D00103',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00103')
        [ig,ig_,_] = s2x_ii('D00203',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00203')
        [ig,ig_,_] = s2x_ii('D00303',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00303')
        [ig,ig_,_] = s2x_ii('D00403',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00403')
        [ig,ig_,_] = s2x_ii('D00104',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00104')
        [ig,ig_,_] = s2x_ii('D00204',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00204')
        [ig,ig_,_] = s2x_ii('D00304',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00304')
        [ig,ig_,_] = s2x_ii('D00404',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00404')
        [ig,ig_,_] = s2x_ii('D00105',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00105')
        [ig,ig_,_] = s2x_ii('D00205',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00205')
        [ig,ig_,_] = s2x_ii('D00305',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00305')
        [ig,ig_,_] = s2x_ii('D00405',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'D00405')
        v_['TM1'] = -1+v_['T']
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            [ig,ig_,_] = s2x_ii('SMOOTH'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SMOOTH'+str(I))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2x_ii('X00101',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00101')
        ig = ig_['K01']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00101']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00201',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00201')
        ig = ig_['K01']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00201']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00301',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00301')
        ig = ig_['K01']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00301']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00401',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00401')
        ig = ig_['K01']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00401']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00102')
        ig = ig_['K02']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00102']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00202',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00202')
        ig = ig_['K02']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00202']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00302',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00302')
        ig = ig_['K02']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00302']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00402',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00402')
        ig = ig_['K02']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00402']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00103',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00103')
        ig = ig_['K03']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00103']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00203')
        ig = ig_['K03']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00203']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00303',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00303')
        ig = ig_['K03']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00303']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00403',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00403')
        ig = ig_['K03']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00403']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00104',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00104')
        ig = ig_['K04']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00104']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00204')
        ig = ig_['K04']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00204']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00304',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00304')
        ig = ig_['K04']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00304']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00404',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00404')
        ig = ig_['K04']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00404']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00105',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00105')
        ig = ig_['K05']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00105']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00205',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00205')
        ig = ig_['K05']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00205']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00305',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00305')
        ig = ig_['K05']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00305']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('X00405',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X00405')
        ig = ig_['K05']
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        ig = ig_['D00405']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00101',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00101')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(1.000000)+pbm.A[ig,iv]
        ig = ig_['D00101']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00101',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00101')
        ig = ig_['D00102']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00101',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00101')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00101']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00201',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00201')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2.000000)+pbm.A[ig,iv]
        ig = ig_['D00201']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00201',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00201')
        ig = ig_['D00202']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00201',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00201')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3)+pbm.A[ig,iv]
        ig = ig_['D00201']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00301',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00301')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3.000000)+pbm.A[ig,iv]
        ig = ig_['D00301']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00301',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00301')
        ig = ig_['D00302']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00301',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00301')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00301']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00401',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00401')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(4.000000)+pbm.A[ig,iv]
        ig = ig_['D00401']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00401',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00401')
        ig = ig_['D00402']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00401',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00401')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(5)+pbm.A[ig,iv]
        ig = ig_['D00401']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00102')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(1.000000)+pbm.A[ig,iv]
        ig = ig_['D00102']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00102')
        ig = ig_['D00103']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00102',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00102')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00102']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00202',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00202')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2.000000)+pbm.A[ig,iv]
        ig = ig_['D00202']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00202',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00202')
        ig = ig_['D00203']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00202',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00202')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3)+pbm.A[ig,iv]
        ig = ig_['D00202']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00302',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00302')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3.000000)+pbm.A[ig,iv]
        ig = ig_['D00302']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00302',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00302')
        ig = ig_['D00303']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00302',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00302')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00302']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00402',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00402')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(4.000000)+pbm.A[ig,iv]
        ig = ig_['D00402']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00402',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00402')
        ig = ig_['D00403']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00402',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00402')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(5)+pbm.A[ig,iv]
        ig = ig_['D00402']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00103',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00103')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(1.000000)+pbm.A[ig,iv]
        ig = ig_['D00103']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00103',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00103')
        ig = ig_['D00104']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00103',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00103')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00103']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00203')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2.000000)+pbm.A[ig,iv]
        ig = ig_['D00203']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00203')
        ig = ig_['D00204']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00203',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00203')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3)+pbm.A[ig,iv]
        ig = ig_['D00203']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00303',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00303')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3.000000)+pbm.A[ig,iv]
        ig = ig_['D00303']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00303',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00303')
        ig = ig_['D00304']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00303',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00303')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00303']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00403',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00403')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(4.000000)+pbm.A[ig,iv]
        ig = ig_['D00403']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00403',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00403')
        ig = ig_['D00404']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00403',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00403')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(5)+pbm.A[ig,iv]
        ig = ig_['D00403']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00104',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00104')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(1.000000)+pbm.A[ig,iv]
        ig = ig_['D00104']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00104',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00104')
        ig = ig_['D00105']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00104',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00104')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00104']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00204')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2.000000)+pbm.A[ig,iv]
        ig = ig_['D00204']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00204')
        ig = ig_['D00205']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00204',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00204')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3)+pbm.A[ig,iv]
        ig = ig_['D00204']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00304',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00304')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3.000000)+pbm.A[ig,iv]
        ig = ig_['D00304']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00304',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00304')
        ig = ig_['D00305']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00304',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00304')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00304']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00404',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00404')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(4.000000)+pbm.A[ig,iv]
        ig = ig_['D00404']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00404',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00404')
        ig = ig_['D00405']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00404',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00404')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(5)+pbm.A[ig,iv]
        ig = ig_['D00404']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00105',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00105')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(1.000000)+pbm.A[ig,iv]
        ig = ig_['D00105']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00105',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00105')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00105']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00205',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00205')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2.000000)+pbm.A[ig,iv]
        ig = ig_['D00205']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00205',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00205')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3)+pbm.A[ig,iv]
        ig = ig_['D00205']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00305',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00305')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(3.000000)+pbm.A[ig,iv]
        ig = ig_['D00305']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00305',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00305')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(2)+pbm.A[ig,iv]
        ig = ig_['D00305']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('I00405',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I00405')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(4.000000)+pbm.A[ig,iv]
        ig = ig_['D00405']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [iv,ix_,_] = s2x_ii('Y00405',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Y00405')
        ig = ig_['COST']
        pbm.A[ig,iv] = float(5)+pbm.A[ig,iv]
        ig = ig_['D00405']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['K01'],float(3))
        pbm.gconst = arrset(pbm.gconst,ig_['K02'],float(6))
        pbm.gconst = arrset(pbm.gconst,ig_['K03'],float(10))
        pbm.gconst = arrset(pbm.gconst,ig_['K04'],float(2000))
        pbm.gconst = arrset(pbm.gconst,ig_['K05'],float(18))
        pbm.gconst = arrset(pbm.gconst,ig_['D00101'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00201'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00301'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00401'],float(1.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00102'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00202'],float(1.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00302'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00402'],float(3.333))
        pbm.gconst = arrset(pbm.gconst,ig_['D00103'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00203'],float(2.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00303'],float(3.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00403'],float(3.000))
        pbm.gconst = arrset(pbm.gconst,ig_['D00104'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00204'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00304'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00404'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00105'],float(2.667))
        pbm.gconst = arrset(pbm.gconst,ig_['D00205'],float(2.333))
        pbm.gconst = arrset(pbm.gconst,ig_['D00305'],float(2.333))
        pbm.gconst = arrset(pbm.gconst,ig_['D00405'],float(2.333))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQMRSQ', iet_)
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
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            v_['IP1'] = 1+I
            ename = 'NLE'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQMRSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQMRSQ"])
            pb.x0 = np.zeros((pb.n,1))
            vname = 'X0010'+str(int(v_['IP1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0020'+str(int(v_['IP1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0030'+str(int(v_['IP1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0040'+str(int(v_['IP1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0010'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0020'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0030'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V7')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0040'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V8')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['TM1'])+1):
            ig = ig_['SMOOTH'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['NLE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "LQR2-RY-60-29"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,2.0)
        pbm.efpar = arrset( pbm.efpar,1,0.1)
        pbm.efpar = arrset( pbm.efpar,2,pbm.efpar[1]*pbm.efpar[1])
        return pbm

    @staticmethod
    def eSQMRSQ(pbm,nargout,*args):

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
        f_   = U1MU2**2-pbm.efpar[2]*IV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.efpar[0]*U1MU2
            g_[1] = -pbm.efpar[0]*U1MU2-pbm.efpar[0]*pbm.efpar[2]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = pbm.efpar[0]
                H_[0,1] = -pbm.efpar[0]
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.efpar[0]*(1.0-pbm.efpar[2])
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

