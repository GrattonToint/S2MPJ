from s2mpjlib import *
class  SSEBNLN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SSEBNLN
#    *********
#    The Power Generation problem for the SSGB.
# 
#    Source:
#    N. Gould, private communication.
# 
#    SIF input: Nick Gould, 23 October 1989
# 
#    classification = "LQR2-RN-194-96"
# 
#    period is the number of time periods
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SSEBNLN'

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
        v_['HOURS'] = 24
        v_['DAYS'] = 1
        v_['PERIOD'] = v_['HOURS']*v_['DAYS']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['Z1'] = 517.0
        v_['D1,1'] = 578.0
        v_['D1,2'] = 517.0
        v_['D1,3'] = 461.0
        v_['D1,4'] = 369.0
        v_['D1,5'] = 299.0
        v_['D1,6'] = 269.0
        v_['D1,7'] = 370.0
        v_['D1,8'] = 559.0
        v_['D1,9'] = 689.0
        v_['D1,10'] = 728.0
        v_['D1,11'] = 683.0
        v_['D1,12'] = 626.0
        v_['D1,13'] = 619.0
        v_['D1,14'] = 586.0
        v_['D1,15'] = 582.0
        v_['D1,16'] = 625.0
        v_['D1,17'] = 821.0
        v_['D1,18'] = 883.0
        v_['D1,19'] = 768.0
        v_['D1,20'] = 711.0
        v_['D1,21'] = 677.0
        v_['D1,22'] = 630.0
        v_['D1,23'] = 545.0
        v_['D1,24'] = 565.0
        v_['Z2'] = 400.0
        v_['D2,1'] = 631.0
        v_['D2,2'] = 574.0
        v_['D2,3'] = 521.0
        v_['D2,4'] = 446.0
        v_['D2,5'] = 359.0
        v_['D2,6'] = 336.0
        v_['D2,7'] = 420.0
        v_['D2,8'] = 588.0
        v_['D2,9'] = 697.0
        v_['D2,10'] = 732.0
        v_['D2,11'] = 713.0
        v_['D2,12'] = 682.0
        v_['D2,13'] = 695.0
        v_['D2,14'] = 651.0
        v_['D2,15'] = 645.0
        v_['D2,16'] = 664.0
        v_['D2,17'] = 816.0
        v_['D2,18'] = 858.0
        v_['D2,19'] = 760.0
        v_['D2,20'] = 700.0
        v_['D2,21'] = 659.0
        v_['D2,22'] = 623.0
        v_['D2,23'] = 517.0
        v_['D2,24'] = 542.0
        v_['Z3'] = 1017.0
        v_['D3,1'] = 582.0
        v_['D3,2'] = 501.0
        v_['D3,3'] = 443.0
        v_['D3,4'] = 367.0
        v_['D3,5'] = 288.0
        v_['D3,6'] = 265.0
        v_['D3,7'] = 349.0
        v_['D3,8'] = 503.0
        v_['D3,9'] = 663.0
        v_['D3,10'] = 651.0
        v_['D3,11'] = 625.0
        v_['D3,12'] = 596.0
        v_['D3,13'] = 608.0
        v_['D3,14'] = 566.0
        v_['D3,15'] = 555.0
        v_['D3,16'] = 584.0
        v_['D3,17'] = 763.0
        v_['D3,18'] = 803.0
        v_['D3,19'] = 710.0
        v_['D3,20'] = 648.0
        v_['D3,21'] = 626.0
        v_['D3,22'] = 590.0
        v_['D3,23'] = 486.0
        v_['D3,24'] = 540.0
        v_['Z4'] = 667.0
        v_['D4,1'] = 602.0
        v_['D4,2'] = 533.0
        v_['D4,3'] = 450.0
        v_['D4,4'] = 378.0
        v_['D4,5'] = 298.0
        v_['D4,6'] = 272.0
        v_['D4,7'] = 369.0
        v_['D4,8'] = 539.0
        v_['D4,9'] = 647.0
        v_['D4,10'] = 652.0
        v_['D4,11'] = 607.0
        v_['D4,12'] = 585.0
        v_['D4,13'] = 587.0
        v_['D4,14'] = 549.0
        v_['D4,15'] = 535.0
        v_['D4,16'] = 564.0
        v_['D4,17'] = 748.0
        v_['D4,18'] = 808.0
        v_['D4,19'] = 710.0
        v_['D4,20'] = 646.0
        v_['D4,21'] = 620.0
        v_['D4,22'] = 581.0
        v_['D4,23'] = 483.0
        v_['D4,24'] = 514.0
        v_['Z5'] = 600.0
        v_['D5,1'] = 579.0
        v_['D5,2'] = 518.0
        v_['D5,3'] = 447.0
        v_['D5,4'] = 355.0
        v_['D5,5'] = 284.0
        v_['D5,6'] = 261.0
        v_['D5,7'] = 348.0
        v_['D5,8'] = 530.0
        v_['D5,9'] = 644.0
        v_['D5,10'] = 648.0
        v_['D5,11'] = 607.0
        v_['D5,12'] = 570.0
        v_['D5,13'] = 577.0
        v_['D5,14'] = 536.0
        v_['D5,15'] = 544.0
        v_['D5,16'] = 554.0
        v_['D5,17'] = 716.0
        v_['D5,18'] = 765.0
        v_['D5,19'] = 676.0
        v_['D5,20'] = 631.0
        v_['D5,21'] = 576.0
        v_['D5,22'] = 528.0
        v_['D5,23'] = 445.0
        v_['D5,24'] = 520.0
        v_['Z6'] = 421.0
        v_['D6,1'] = 618.0
        v_['D6,2'] = 547.0
        v_['D6,3'] = 430.0
        v_['D6,4'] = 327.0
        v_['D6,5'] = 249.0
        v_['D6,6'] = 211.0
        v_['D6,7'] = 227.0
        v_['D6,8'] = 258.0
        v_['D6,9'] = 347.0
        v_['D6,10'] = 491.0
        v_['D6,11'] = 524.0
        v_['D6,12'] = 492.0
        v_['D6,13'] = 467.0
        v_['D6,14'] = 418.0
        v_['D6,15'] = 358.0
        v_['D6,16'] = 378.0
        v_['D6,17'] = 544.0
        v_['D6,18'] = 666.0
        v_['D6,19'] = 589.0
        v_['D6,20'] = 533.0
        v_['D6,21'] = 494.0
        v_['D6,22'] = 460.0
        v_['D6,23'] = 404.0
        v_['D6,24'] = 512.0
        v_['Z7'] = 425.0
        v_['D7,1'] = 615.0
        v_['D7,2'] = 587.0
        v_['D7,3'] = 450.0
        v_['D7,4'] = 320.0
        v_['D7,5'] = 235.0
        v_['D7,6'] = 198.0
        v_['D7,7'] = 195.0
        v_['D7,8'] = 173.0
        v_['D7,9'] = 197.0
        v_['D7,10'] = 349.0
        v_['D7,11'] = 441.0
        v_['D7,12'] = 459.0
        v_['D7,13'] = 485.0
        v_['D7,14'] = 445.0
        v_['D7,15'] = 410.0
        v_['D7,16'] = 421.0
        v_['D7,17'] = 568.0
        v_['D7,18'] = 643.0
        v_['D7,19'] = 596.0
        v_['D7,20'] = 566.0
        v_['D7,21'] = 541.0
        v_['D7,22'] = 532.0
        v_['D7,23'] = 454.0
        v_['D7,24'] = 511.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('V'+str(int(v_['0']))+','+str(int(v_['HOURS'])),ix_)
        pb.xnames = (
             arrset(pb.xnames,iv,'V'+str(int(v_['0']))+','+str(int(v_['HOURS']))))
        [iv,ix_,_] = s2mpj_ii('R'+str(int(v_['0']))+','+str(int(v_['HOURS'])),ix_)
        pb.xnames = (
             arrset(pb.xnames,iv,'R'+str(int(v_['0']))+','+str(int(v_['HOURS']))))
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                [iv,ix_,_] = s2mpj_ii('P1'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'P1'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('P2'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'P2'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('QH'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'QH'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('S'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'S'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('QG'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'QG'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('QP'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'QP'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('V'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'V'+str(ID)+','+str(IH))
                [iv,ix_,_] = s2mpj_ii('R'+str(ID)+','+str(IH),ix_)
                pb.xnames=arrset(pb.xnames,iv,'R'+str(ID)+','+str(IH))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['P1'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1000.0)+pbm.A[ig,iv]
                iv = ix_['P2'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1500.0)+pbm.A[ig,iv]
                iv = ix_['QH'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1200.0)+pbm.A[ig,iv]
                iv = ix_['S'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1200.0)+pbm.A[ig,iv]
                iv = ix_['QG'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1200.0)+pbm.A[ig,iv]
                iv = ix_['QP'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(-1200.0)+pbm.A[ig,iv]
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            v_['P'] = -1+ID
            [ig,ig_,_] = s2mpj_ii('H'+str(ID)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'H'+str(ID)+','+str(int(v_['1'])))
            iv = ix_['V'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['V'+str(int(v_['P']))+','+str(int(v_['HOURS']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('H'+str(ID)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'H'+str(ID)+','+str(int(v_['1'])))
            iv = ix_['S'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['QH'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['2']),int(v_['HOURS'])+1):
                v_['IH-1'] = -1+IH
                [ig,ig_,_] = s2mpj_ii('H'+str(ID)+','+str(IH),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'H'+str(ID)+','+str(IH))
                iv = ix_['V'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['V'+str(ID)+','+str(int(v_['IH-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['S'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['QH'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            v_['P'] = -1+ID
            [ig,ig_,_] = s2mpj_ii('R'+str(ID)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(ID)+','+str(int(v_['1'])))
            iv = ix_['R'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['R'+str(int(v_['P']))+','+str(int(v_['HOURS']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('R'+str(ID)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(ID)+','+str(int(v_['1'])))
            iv = ix_['QG'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['QP'+str(ID)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['2']),int(v_['HOURS'])+1):
                v_['IH-1'] = -1+IH
                [ig,ig_,_] = s2mpj_ii('R'+str(ID)+','+str(IH),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'R'+str(ID)+','+str(IH))
                iv = ix_['R'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['R'+str(ID)+','+str(int(v_['IH-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['QG'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['QP'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                [ig,ig_,_] = s2mpj_ii('D'+str(ID)+','+str(IH),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'D'+str(ID)+','+str(IH))
                iv = ix_['P1'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['P2'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['QH'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['QG'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['QP'+str(ID)+','+str(IH)]
                pbm.A[ig,iv] = float(-1.33)+pbm.A[ig,iv]
        for D in range(int(v_['1']),int(v_['DAYS'])+1):
            for H in range(int(v_['1']),int(v_['HOURS'])+1):
                [ig,ig_,_] = s2mpj_ii('QG*QP'+str(D)+','+str(H),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'QG*QP'+str(D)+','+str(H))
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
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['H'+str(ID)+','+str(IH)],float(v_['Z'+str(ID)])))
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            v_['0.01Z'] = 0.01*v_['Z'+str(ID)]
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['R'+str(ID)+','+str(IH)],float(v_['0.01Z'])))
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['D'+str(ID)+','+str(IH)],float(v_['D'+str(ID)+','+str(IH)])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['V'+str(int(v_['0']))+','+str(int(v_['HOURS']))]] = 240000.0
        pb.xupper[ix_['V'+str(int(v_['0']))+','+str(int(v_['HOURS']))]] = 240000.0
        pb.xlower[ix_['R'+str(int(v_['0']))+','+str(int(v_['HOURS']))]] = 3500.0
        pb.xupper[ix_['R'+str(int(v_['0']))+','+str(int(v_['HOURS']))]] = 3500.0
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                pb.xlower[ix_['P1'+str(ID)+','+str(IH)]] = 70.0
                pb.xlower[ix_['P2'+str(ID)+','+str(IH)]] = 90.0
                pb.xlower[ix_['QH'+str(ID)+','+str(IH)]] = 25.0
                pb.xlower[ix_['V'+str(ID)+','+str(IH)]] = 180000.0
                pb.xupper[ix_['P1'+str(ID)+','+str(IH)]] = 325.0
                pb.xupper[ix_['P2'+str(ID)+','+str(IH)]] = 290.0
                pb.xupper[ix_['QH'+str(ID)+','+str(IH)]] = 500.0
                pb.xupper[ix_['QP'+str(ID)+','+str(IH)]] = 225.0
                pb.xupper[ix_['QG'+str(ID)+','+str(IH)]] = 300.0
                pb.xupper[ix_['V'+str(ID)+','+str(IH)]] = 280000.0
                pb.xupper[ix_['R'+str(ID)+','+str(IH)]] = 6000.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('V'+str(int(v_['0']))+','+str(int(v_['HOURS'])) in ix_):
            pb.x0[ix_['V'+str(int(v_['0']))+','+str(int(v_['HOURS']))]]  = (
                  float(240000.0))
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['V'+str(int(v_['0']))+','+str(int(v_['HOURS']))]),float(240000.0)))
        if('R'+str(int(v_['0']))+','+str(int(v_['HOURS'])) in ix_):
            pb.x0[ix_['R'+str(int(v_['0']))+','+str(int(v_['HOURS']))]] = float(3500.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['R'+str(int(v_['0']))+','+str(int(v_['HOURS']))]),float(3500.0)))
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                if('P1'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['P1'+str(ID)+','+str(IH)]] = float(70.0)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1'+str(ID)+','+str(IH)]),float(70.0)))
                if('P2'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['P2'+str(ID)+','+str(IH)]] = float(90.0)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2'+str(ID)+','+str(IH)]),float(90.0)))
                if('QH'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['QH'+str(ID)+','+str(IH)]] = float(25.0)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['QH'+str(ID)+','+str(IH)]),float(25.0)))
                if('QP'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['QP'+str(ID)+','+str(IH)]] = float(225.0)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['QP'+str(ID)+','+str(IH)]),float(225.0)))
                if('V'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['V'+str(ID)+','+str(IH)]] = float(240000.0)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['V'+str(ID)+','+str(IH)]),float(240000.0)))
                if('R'+str(ID)+','+str(IH) in ix_):
                    pb.x0[ix_['R'+str(ID)+','+str(IH)]] = float(3500)
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['R'+str(ID)+','+str(IH)]),float(3500)))
        pass
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'QP')
        elftv = loaset(elftv,it,1,'QG')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for ID in range(int(v_['1']),int(v_['DAYS'])+1):
            for IH in range(int(v_['1']),int(v_['HOURS'])+1):
                ename = 'P'+str(ID)+','+str(IH)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                vname = 'QP'+str(ID)+','+str(IH)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='QP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'QG'+str(ID)+','+str(IH)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='QG')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for D in range(int(v_['1']),int(v_['DAYS'])+1):
            for H in range(int(v_['1']),int(v_['HOURS'])+1):
                ig = ig_['QG*QP'+str(D)+','+str(H)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(D)+','+str(H)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1.617060D+07
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-RN-194-96"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[1]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[1] = EV_[0]
            g_[0] = EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = 1.0e+0
                H_[0,1] = H_[1,0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

