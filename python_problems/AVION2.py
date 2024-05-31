from s2xlib import *
class  AVION2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : AVION2
#    *********
# 
#    Dassault France avion (airplane design) problem
# 
#    SIF input:  A. R. Conn, June 1993.
# 
#    classification = "OLR2-RN-49-15"
# 
#    Define useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'AVION2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'AVION2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['790'] = 790.0
        v_['1/790'] = 1.0/v_['790']
        v_['24000'] = 24000.0
        v_['1/24000'] = 1.0/v_['24000']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('SR',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SR')
        [iv,ix_,_] = s2x_ii('LR',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LR')
        [iv,ix_,_] = s2x_ii('PK',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PK')
        [iv,ix_,_] = s2x_ii('EF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EF')
        [iv,ix_,_] = s2x_ii('SX',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SX')
        [iv,ix_,_] = s2x_ii('LX',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LX')
        [iv,ix_,_] = s2x_ii('SD',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SD')
        [iv,ix_,_] = s2x_ii('SK',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SK')
        [iv,ix_,_] = s2x_ii('ST',ix_)
        pb.xnames=arrset(pb.xnames,iv,'ST')
        [iv,ix_,_] = s2x_ii('SF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SF')
        [iv,ix_,_] = s2x_ii('LF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'LF')
        [iv,ix_,_] = s2x_ii('AM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'AM')
        [iv,ix_,_] = s2x_ii('CA',ix_)
        pb.xnames=arrset(pb.xnames,iv,'CA')
        [iv,ix_,_] = s2x_ii('CB',ix_)
        pb.xnames=arrset(pb.xnames,iv,'CB')
        [iv,ix_,_] = s2x_ii('SO',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SO')
        [iv,ix_,_] = s2x_ii('SS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'SS')
        [iv,ix_,_] = s2x_ii('IMPDER',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPDER')
        [iv,ix_,_] = s2x_ii('IMPK',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPK')
        [iv,ix_,_] = s2x_ii('IMPFUS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPFUS')
        [iv,ix_,_] = s2x_ii('QI',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QI')
        [iv,ix_,_] = s2x_ii('PT',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PT')
        [iv,ix_,_] = s2x_ii('MV',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MV')
        [iv,ix_,_] = s2x_ii('MC',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MC')
        [iv,ix_,_] = s2x_ii('MD',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MD')
        [iv,ix_,_] = s2x_ii('PD',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PD')
        [iv,ix_,_] = s2x_ii('NS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NS')
        [iv,ix_,_] = s2x_ii('VS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VS')
        [iv,ix_,_] = s2x_ii('CR',ix_)
        pb.xnames=arrset(pb.xnames,iv,'CR')
        [iv,ix_,_] = s2x_ii('PM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'PM')
        [iv,ix_,_] = s2x_ii('DV',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DV')
        [iv,ix_,_] = s2x_ii('MZ',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MZ')
        [iv,ix_,_] = s2x_ii('VN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'VN')
        [iv,ix_,_] = s2x_ii('QV',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QV')
        [iv,ix_,_] = s2x_ii('QF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'QF')
        [iv,ix_,_] = s2x_ii('IMPTRAIN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPTRAIN')
        [iv,ix_,_] = s2x_ii('IMPMOT',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPMOT')
        [iv,ix_,_] = s2x_ii('IMPNMOT',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPNMOT')
        [iv,ix_,_] = s2x_ii('IMPPET',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPPET')
        [iv,ix_,_] = s2x_ii('IMPPIL',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPPIL')
        [iv,ix_,_] = s2x_ii('IMPCAN',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPCAN')
        [iv,ix_,_] = s2x_ii('IMPSNA',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPSNA')
        [iv,ix_,_] = s2x_ii('MS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MS')
        [iv,ix_,_] = s2x_ii('EL',ix_)
        pb.xnames=arrset(pb.xnames,iv,'EL')
        [iv,ix_,_] = s2x_ii('DE',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DE')
        [iv,ix_,_] = s2x_ii('DS',ix_)
        pb.xnames=arrset(pb.xnames,iv,'DS')
        [iv,ix_,_] = s2x_ii('IMPVOIL',ix_)
        pb.xnames=arrset(pb.xnames,iv,'IMPVOIL')
        [iv,ix_,_] = s2x_ii('NM',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NM')
        [iv,ix_,_] = s2x_ii('NP',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NP')
        [iv,ix_,_] = s2x_ii('NG',ix_)
        pb.xnames=arrset(pb.xnames,iv,'NG')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('E1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E1')
        iv = ix_['SD']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SR']
        pbm.A[ig,iv] = float(-0.13)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E2')
        iv = ix_['SX']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SR']
        pbm.A[ig,iv] = float(-0.7)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E3')
        iv = ix_['LX']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['LR']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['SK']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E5',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E5')
        iv = ix_['SF']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['ST']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['SD']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['SX']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['SK']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E6',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['CA']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E7',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['AM']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        iv = ix_['SO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SS']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E8',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['AM']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E9',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['IMPDER']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SD']
        pbm.A[ig,iv] = float(-27.5)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E10',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['IMPK']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SK']
        pbm.A[ig,iv] = float(-70.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E11',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E11')
        iv = ix_['IMPFUS']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['SF']
        pbm.A[ig,iv] = float(-20.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E12',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E12')
        iv = ix_['MD']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['MV']
        pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E13',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['QI']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E14',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['PT']
        pbm.A[ig,iv] = float(1000.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E15',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E15')
        iv = ix_['QF']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['QI']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['QV']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E16',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['VN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['VS']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['QF']
        pbm.A[ig,iv] = float(v_['1/790'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E17',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E17')
        iv = ix_['IMPTRAIN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['MV']
        pbm.A[ig,iv] = float(-0.137)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E18',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['IMPMOT']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E19',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E19')
        iv = ix_['IMPNMOT']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NM']
        pbm.A[ig,iv] = float(-35.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E20',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E20')
        iv = ix_['IMPPET']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['QI']
        pbm.A[ig,iv] = float(-0.043)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E21',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E21')
        iv = ix_['IMPPIL']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NP']
        pbm.A[ig,iv] = float(-200.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E22',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E22')
        iv = ix_['IMPCAN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NG']
        pbm.A[ig,iv] = float(-120.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E23',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E23')
        iv = ix_['IMPSNA']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['NS']
        pbm.A[ig,iv] = float(-300.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E24',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E24')
        iv = ix_['MC']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['MV']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['NP']
        pbm.A[ig,iv] = float(95.0)+pbm.A[ig,iv]
        iv = ix_['NG']
        pbm.A[ig,iv] = float(70.0)+pbm.A[ig,iv]
        iv = ix_['NM']
        pbm.A[ig,iv] = float(660.0)+pbm.A[ig,iv]
        iv = ix_['QI']
        pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E25',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'E25')
        iv = ix_['MZ']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IMPTRAIN']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['IMPNMOT']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IMPPET']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IMPPIL']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IMPCAN']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['IMPSNA']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E26',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['ST']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E27',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['SR']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E28',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['QV']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E29',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['SO']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E30',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['SS']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E31',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['CB']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E32',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['IMPVOIL']
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
        pbm.gconst = arrset(pbm.gconst,ig_['E13'],float(1000.0))
        pbm.gconst = arrset(pbm.gconst,ig_['E16'],float(-2.0))
        pbm.gconst = arrset(pbm.gconst,ig_['E23'],float(400.0))
        pbm.gconst = arrset(pbm.gconst,ig_['E24'],float(380.0))
        pbm.gconst = arrset(pbm.gconst,ig_['E25'],float(-290.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['SR']] = 10.0
        pb.xlower[ix_['LR']] = 0.0
        pb.xlower[ix_['PK']] = 0.0
        pb.xlower[ix_['EF']] = 0.0
        pb.xlower[ix_['SX']] = 7.0
        pb.xlower[ix_['LX']] = 1.5
        pb.xlower[ix_['SD']] = 2.0
        pb.xlower[ix_['SK']] = 2.0
        pb.xlower[ix_['ST']] = 30.0
        pb.xlower[ix_['SF']] = 20.0
        pb.xlower[ix_['LF']] = 0.001
        pb.xlower[ix_['AM']] = 0.0
        pb.xlower[ix_['CA']] = -0.2
        pb.xlower[ix_['CB']] = 0.1
        pb.xlower[ix_['SO']] = 0.0
        pb.xlower[ix_['SS']] = 0.0
        pb.xlower[ix_['IMPDER']] = 100.0
        pb.xlower[ix_['IMPK']] = 500.0
        pb.xlower[ix_['IMPFUS']] = 500.0
        pb.xlower[ix_['QI']] = 1000.0
        pb.xlower[ix_['PT']] = 2.0
        pb.xlower[ix_['MV']] = 2000.0
        pb.xlower[ix_['MC']] = 3000.0
        pb.xlower[ix_['MD']] = 5000.0
        pb.xlower[ix_['PD']] = 0.2
        pb.xlower[ix_['NS']] = 1.0
        pb.xlower[ix_['VS']] = 0.0
        pb.xlower[ix_['CR']] = 100.0
        pb.xlower[ix_['PM']] = 4.0
        pb.xlower[ix_['DV']] = 0.0
        pb.xlower[ix_['MZ']] = 500.0
        pb.xlower[ix_['VN']] = 10.0
        pb.xlower[ix_['QV']] = 250.0
        pb.xlower[ix_['QF']] = 750.0
        pb.xlower[ix_['IMPTRAIN']] = 250.0
        pb.xlower[ix_['IMPMOT']] = 10.0
        pb.xlower[ix_['IMPNMOT']] = 35.0
        pb.xlower[ix_['IMPPET']] = 100.0
        pb.xlower[ix_['IMPPIL']] = 200.0
        pb.xlower[ix_['IMPCAN']] = 120.0
        pb.xlower[ix_['IMPSNA']] = 700.0
        pb.xlower[ix_['MS']] = 100.0
        pb.xlower[ix_['EL']] = 2.0
        pb.xlower[ix_['DE']] = 0.0
        pb.xlower[ix_['DS']] = 0.0
        pb.xlower[ix_['IMPVOIL']] = 500.0
        pb.xupper[ix_['SR']] = 150.0
        pb.xupper[ix_['LR']] = 10.0
        pb.xupper[ix_['PK']] = 10.0
        pb.xupper[ix_['EF']] = 5.0
        pb.xupper[ix_['SX']] = 120.0
        pb.xupper[ix_['LX']] = 8.0
        pb.xupper[ix_['SD']] = 20.0
        pb.xupper[ix_['SK']] = 30.0
        pb.xupper[ix_['ST']] = 500.0
        pb.xupper[ix_['SF']] = 200.0
        pb.xupper[ix_['LF']] = 20.0
        pb.xupper[ix_['AM']] = 10.0
        pb.xupper[ix_['CA']] = -0.001
        pb.xupper[ix_['CB']] = 2.0
        pb.xupper[ix_['SO']] = 1.0
        pb.xupper[ix_['SS']] = 2.0
        pb.xupper[ix_['IMPDER']] = 1000.0
        pb.xupper[ix_['IMPK']] = 5000.0
        pb.xupper[ix_['IMPFUS']] = 5000.0
        pb.xupper[ix_['QI']] = 20000.0
        pb.xupper[ix_['PT']] = 30.0
        pb.xupper[ix_['MV']] = 20000.0
        pb.xupper[ix_['MC']] = 30000.0
        pb.xupper[ix_['MD']] = 50000.0
        pb.xupper[ix_['PD']] = 0.8
        pb.xupper[ix_['NS']] = 5.0
        pb.xupper[ix_['VS']] = 20.0
        pb.xupper[ix_['CR']] = 400.0
        pb.xupper[ix_['PM']] = 15.0
        pb.xupper[ix_['DV']] = 10.0
        pb.xupper[ix_['MZ']] = 10000.0
        pb.xupper[ix_['VN']] = 50.0
        pb.xupper[ix_['QV']] = 5000.0
        pb.xupper[ix_['QF']] = 15000.0
        pb.xupper[ix_['IMPTRAIN']] = 3000.0
        pb.xupper[ix_['IMPMOT']] = 5000.0
        pb.xupper[ix_['IMPNMOT']] = 70.0
        pb.xupper[ix_['IMPPET']] = 3000.0
        pb.xupper[ix_['IMPPIL']] = 400.0
        pb.xupper[ix_['IMPCAN']] = 240.0
        pb.xupper[ix_['IMPSNA']] = 1900.0
        pb.xupper[ix_['MS']] = 1000.0
        pb.xupper[ix_['EL']] = 20.0
        pb.xupper[ix_['DE']] = 1.0
        pb.xupper[ix_['DS']] = 2.0
        pb.xupper[ix_['IMPVOIL']] = 5000.0
        pb.xupper[ix_['NM']] = 2.0
        pb.xlower[ix_['NM']] = 1.0
        pb.xupper[ix_['NP']] = 2.0
        pb.xlower[ix_['NP']] = 1.0
        pb.xupper[ix_['NG']] = 2.0
        pb.xlower[ix_['NG']] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['SR']] = float(2.7452e+01)
        pb.x0[ix_['LR']] = float(1.5000)
        pb.x0[ix_['PK']] = float(1.0000e+01)
        pb.x0[ix_['EF']] = float(0.0000)
        pb.x0[ix_['SX']] = float(1.9217e+01)
        pb.x0[ix_['LX']] = float(1.5000)
        pb.x0[ix_['SD']] = float(3.5688)
        pb.x0[ix_['SK']] = float(4.0696)
        pb.x0[ix_['ST']] = float(3.4315e+01)
        pb.x0[ix_['SF']] = float(8.8025e+01)
        pb.x0[ix_['LF']] = float(5.1306)
        pb.x0[ix_['AM']] = float(0.0000)
        pb.x0[ix_['CA']] = float(-1.4809e-01)
        pb.x0[ix_['CB']] = float(7.5980e-01)
        pb.x0[ix_['SO']] = float(0.0000)
        pb.x0[ix_['SS']] = float(0.0000)
        pb.x0[ix_['IMPDER']] = float(1.1470e+02)
        pb.x0[ix_['IMPK']] = float(5.0000e+02)
        pb.x0[ix_['IMPFUS']] = float(1.7605e+03)
        pb.x0[ix_['QI']] = float(2.3256e+03)
        pb.x0[ix_['PT']] = float(5.6788)
        pb.x0[ix_['MV']] = float(1.4197e+04)
        pb.x0[ix_['MC']] = float(1.2589e+04)
        pb.x0[ix_['MD']] = float(2.8394e+04)
        pb.x0[ix_['PD']] = float(2.0000e-01)
        pb.x0[ix_['NS']] = float(1.0000)
        pb.x0[ix_['VS']] = float(0.0000)
        pb.x0[ix_['CR']] = float(1.0000e+02)
        pb.x0[ix_['PM']] = float(1.5000e+01)
        pb.x0[ix_['DV']] = float(0.0000)
        pb.x0[ix_['MZ']] = float(5.0000e+02)
        pb.x0[ix_['VN']] = float(1.0000e+01)
        pb.x0[ix_['QV']] = float(8.1490e+02)
        pb.x0[ix_['QF']] = float(3.1405e+03)
        pb.x0[ix_['IMPTRAIN']] = float(1.9450e+03)
        pb.x0[ix_['IMPMOT']] = float(1.9085e+02)
        pb.x0[ix_['IMPNMOT']] = float(3.5000e+01)
        pb.x0[ix_['IMPPET']] = float(1.0000e+02)
        pb.x0[ix_['IMPPIL']] = float(2.0000e+02)
        pb.x0[ix_['IMPCAN']] = float(1.2000e+02)
        pb.x0[ix_['IMPSNA']] = float(7.0000e+02)
        pb.x0[ix_['MS']] = float(1.0000e+03)
        pb.x0[ix_['EL']] = float(4.9367)
        pb.x0[ix_['DE']] = float(0.0000)
        pb.x0[ix_['DS']] = float(0.0000)
        pb.x0[ix_['IMPVOIL']] = float(5.0000e+03)
        pb.x0[ix_['NM']] = float(1.0000)
        pb.x0[ix_['NP']] = float(1.0000)
        pb.x0[ix_['NG']] = float(1.0000)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2x_ii( 'eQDdSQ', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        [it,iet_,_] = s2x_ii( 'en12', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2x_ii( 'en12d1', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2x_ii( 'eQT', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2x_ii( 'en1dLIN', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2x_ii( 'eSQRT', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2x_ii( 'eSURD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2x_ii( 'eSQPRD', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2x_ii( 'eCBdSQQD', iet_)
        elftv = loaset(elftv,it,0,'W')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        elftv = loaset(elftv,it,3,'Z')
        [it,iet_,_] = s2x_ii( 'eSREL', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        elftv = loaset(elftv,it,2,'X')
        elftv = loaset(elftv,it,3,'Y')
        elftv = loaset(elftv,it,4,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'EL1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'PK'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SR'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eQDdSQ')
        ielftype = arrset(ielftype, ie, iet_["eQDdSQ"])
        vname = 'SS'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SO'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'CB'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'LF'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en12')
        ielftype = arrset(ielftype, ie, iet_["en12"])
        vname = 'EF'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'LF'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL4'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en12d1')
        ielftype = arrset(ielftype, ie, iet_["en12d1"])
        vname = 'SO'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'CB'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'CA'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL5'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'SD'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL6'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'SK'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL7'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'MV'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL8'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'MD'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PD'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL9'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eQT')
        ielftype = arrset(ielftype, ie, iet_["eQT"])
        vname = 'MZ'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'CR'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL10'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'DV'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PT'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL11'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en1dLIN')
        ielftype = arrset(ielftype, ie, iet_["en1dLIN"])
        vname = 'PT'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PM'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL12'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQRT')
        ielftype = arrset(ielftype, ie, iet_["eSQRT"])
        vname = 'PT'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL13'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
        ielftype = arrset(ielftype, ie, iet_["en2PR"])
        vname = 'SR'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'NM'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL14'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eQT')
        ielftype = arrset(ielftype, ie, iet_["eQT"])
        vname = 'MD'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MS'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL15'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSURD')
        ielftype = arrset(ielftype, ie, iet_["eSURD"])
        vname = 'SX'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'EL'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'LX'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL16'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQPRD')
        ielftype = arrset(ielftype, ie, iet_["eSQPRD"])
        vname = 'DE'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PT'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL17'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQPRD')
        ielftype = arrset(ielftype, ie, iet_["eSQPRD"])
        vname = 'DS'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'PT'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL18'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eCBdSQQD')
        ielftype = arrset(ielftype, ie, iet_["eCBdSQQD"])
        vname = 'VN'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'CA'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'LF'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SO'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'EL19'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSREL')
        ielftype = arrset(ielftype, ie, iet_["eSREL"])
        vname = 'SX'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'MC'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'LX'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'SR'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'EL'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['E4']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.01))
        ig = ig_['E6']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['E7']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(0.01))
        ig = ig_['E8']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.25))
        ig = ig_['E9']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.3))
        ig = ig_['E10']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(8.6))
        ig = ig_['E13']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/24000']))
        ig = ig_['E14']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['E16']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL9'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL10'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        ig = ig_['E18']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1000.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL12'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-12.0))
        ig = ig_['E26']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL13'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.25))
        ig = ig_['E27']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL14'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['E28']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL15'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.4))
        ig = ig_['E29']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL16'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.785))
        ig = ig_['E30']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL17'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-0.785))
        ig = ig_['E31']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL18'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.0))
        ig = ig_['E32']
        pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL19'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.15))
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
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-RN-49-15"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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
                H_[0,1] = 1.0e0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQDdSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        QD = EV_[0]-EV_[1]-EV_[2]*EV_[3]
        SQ = EV_[3]**2
        RSQ = 1.0e0/SQ
        QDOSQ = QD/SQ
        f_   = QDOSQ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = RSQ
            g_[1] = -RSQ
            g_[2] = -1.0e0/EV_[3]
            g_[3] = -EV_[2]/SQ-2.0e0*QDOSQ/EV_[3]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,3] = -2.0e0/(SQ*EV_[3])
                H_[3,0] = H_[0,3]
                H_[1,3] = 2.0e0/(SQ*EV_[3])
                H_[3,1] = H_[1,3]
                H_[2,3] = RSQ
                H_[3,2] = H_[2,3]
                H_[3,3] = (4.0e0*EV_[2])/(SQ*EV_[3])+6.0e0*QDOSQ/SQ
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en12(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]/EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EV_[1]
            g_[1] = -EV_[0]/EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0e0/EV_[1]**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (2.0e0*EV_[0])/EV_[1]**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en12d1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ZSQ = EV_[2]**2
        YSQ = EV_[1]**2
        f_   = (EV_[0]*YSQ)/EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = YSQ/EV_[2]
            g_[1] = (2.0e0*EV_[0]*EV_[1])/EV_[2]
            g_[2] = -(EV_[0]*YSQ)/ZSQ
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = (2.0e0*EV_[1])/EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = -YSQ/ZSQ
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0e0*EV_[0])/EV_[2]
                H_[1,2] = -(2.0e0*EV_[0]*EV_[1])/ZSQ
                H_[2,1] = H_[1,2]
                H_[2,2] = (2.0e0*EV_[0]*YSQ)/(ZSQ*EV_[2])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]/EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EV_[1]
            g_[1] = -EV_[0]/EV_[1]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0e0/EV_[1]**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (2.0e0*EV_[0])/EV_[1]**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en1dLIN(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LIN = EV_[1]+20.0e0
        f_   = EV_[0]/LIN
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/LIN
            g_[1] = -EV_[0]/LIN**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0e0/LIN**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (2.0e0*EV_[0])/LIN**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQRT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        RTZ = np.sqrt(EV_[0])
        f_   = RTZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.5e0/RTZ
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -0.25e0/(EV_[0]*RTZ)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSURD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        RTX = np.sqrt(EV_[0])
        RTZ = np.sqrt(EV_[2])
        XRTX = EV_[0]*np.sqrt(EV_[0])
        ZRTZ = EV_[2]*np.sqrt(EV_[2])
        f_   = XRTX*EV_[1]/RTZ
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.5e0*RTX*EV_[1]/RTZ
            g_[1] = XRTX/RTZ
            g_[2] = -(0.5e0*XRTX*EV_[1])/ZRTZ
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = 0.75e0*EV_[1]/(RTZ*RTX)
                H_[0,1] = 1.5e0*RTX/RTZ
                H_[1,0] = H_[0,1]
                H_[0,2] = -(0.75e0*RTX*EV_[1])/ZRTZ
                H_[2,0] = H_[0,2]
                H_[1,2] = -(0.5e0*XRTX)/ZRTZ
                H_[2,1] = H_[1,2]
                H_[2,2] = (0.75e0*XRTX*EV_[1])/(ZRTZ*EV_[2])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQPRD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**2*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0e0*EV_[0]*EV_[1]
            g_[1] = EV_[0]**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0e0*EV_[1]
                H_[0,1] = 2.0e0*EV_[0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eCBdSQQD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        YCB = EV_[2]**3
        CB = EV_[0]-EV_[1]*YCB
        TMZY = 3.0e0-EV_[3]*EV_[2]
        SQQD = EV_[2]**2*TMZY
        DYSQQD = 3.0e0*EV_[2]*(2.0e0-EV_[3]*EV_[2])
        f_   = CB/SQQD
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0e0/SQQD
            g_[1] = -YCB/SQQD
            g_[2] = -3.0e0*EV_[1]*EV_[2]**2/SQQD-(CB*DYSQQD)/(SQQD*SQQD)
            g_[3] = (YCB*CB)/SQQD**2
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,2] = -DYSQQD/SQQD**2
                H_[2,0] = H_[0,2]
                H_[0,3] = YCB/SQQD**2
                H_[3,0] = H_[0,3]
                H_[1,2] = -3.0e0*(1.0e0-(TMZY-1.0e0)/TMZY)/TMZY
                H_[2,1] = H_[1,2]
                H_[1,3] = -(YCB**2*CB)/SQQD**2
                H_[3,1] = H_[1,3]
                H_[2,2] = (-6.0e0*EV_[1]*EV_[2]/SQQD+(3.0e0*EV_[1]*EV_[2]**2*DYSQQD)/
                     SQQD**2+(2.0e0*CB*DYSQQD**2)/SQQD**3-(6.0e0*CB*(1.0-EV_[3]*EV_[2]))/SQQD**2)
                H_[2,3] = (-(3.0e0*EV_[1]*EV_[2])/TMZY**2-(6.0e0*EV_[2]**4*(TMZY-1.0e0)*CB)/
                     SQQD**3+3.0e0*CB*EV_[2]**2/SQQD**2)
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0e0*(YCB**2*CB)/SQQD**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSREL(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        SRLIN = 15.0e0+0.15e0*EV_[0]
        SRPD = EV_[1]*EV_[2]
        SRQD = 50.0e0*EV_[3]*EV_[4]
        SRQT = SRPD/SRQD
        SRRT = np.sqrt(SRQT)
        SRSUM = 15.0e0+0.3e0*EV_[0]
        SRVLN = EV_[0]*SRLIN
        f_   = EV_[0]*SRLIN*(SRQT*SRRT+8.0e0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SRSUM*(SRQT*SRRT+8.0e0)
            g_[1] = 1.5e0*(SRVLN*SRRT*EV_[2]/SRQD)
            g_[2] = 1.5e0*(SRVLN*SRRT*EV_[1]/SRQD)
            g_[3] = -1.5e0*SRVLN*SRRT*SRQT/EV_[3]
            g_[4] = -1.5e0*SRVLN*SRRT*SRQT/EV_[4]
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = 0.3e0*(SRQT*SRRT+8.0)
                H_[0,1] = 1.5e0*(SRSUM*SRRT*EV_[2]/SRQD)
                H_[1,0] = H_[0,1]
                H_[0,2] = 1.5e0*(SRSUM*SRRT*EV_[1]/SRQD)
                H_[2,0] = H_[0,2]
                H_[0,3] = -1.5e0*SRSUM*SRRT*SRQT/EV_[3]
                H_[3,0] = H_[0,3]
                H_[0,4] = -1.5e0*SRSUM*SRRT*SRQT/EV_[4]
                H_[4,0] = H_[0,4]
                H_[1,1] = (0.75e0*SRVLN*EV_[2]**2)/(SRQD**2*SRRT)
                H_[1,2] = SRVLN*((0.75e0*SRPD)/(SRQD**2*SRRT)+(1.5e0*SRRT)/SRQD)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[1]*EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[1]*EV_[4])
                H_[4,1] = H_[1,4]
                H_[2,2] = (SRVLN*0.75e0*EV_[1]*EV_[1])/(SRRT*SRQD**2)
                H_[2,3] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[2]*EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = -(SRVLN*2.25e0*SRRT*SRQT)/(EV_[2]*EV_[4])
                H_[4,2] = H_[2,4]
                H_[3,3] = (SRVLN*3.75e0*SRRT*SRQT)/EV_[3]**2
                H_[3,4] = (SRVLN*2.25e0*SRRT*SRQT)/(EV_[3]*EV_[4])
                H_[4,3] = H_[3,4]
                H_[4,4] = (SRVLN*3.75e0*SRRT*SRQT)/EV_[4]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQUARE(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

