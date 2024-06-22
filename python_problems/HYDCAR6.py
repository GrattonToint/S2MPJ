from s2mpjlib import *
class  HYDCAR6(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HYDCAR6
#    *********
# 
#    The hydrocarbon-6 problem by Fletcher.
# 
#    Source: Problem 2a in
#    J.J. More',"A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer Seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input : N. Gould, Dec 1989
# 
#    classification = "NOR2-AN-29-29"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HYDCAR6'

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
        v_['0'] = 0
        v_['1'] = 1
        v_['N'] = 6
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        v_['M'] = 3
        v_['K'] = 2
        v_['K-'] = -1+v_['K']
        v_['K+'] = 1+v_['K']
        v_['A1'] = 9.647
        v_['B1'] = -2998.00
        v_['C1'] = 230.66
        v_['A2'] = 9.953
        v_['B2'] = -3448.10
        v_['C2'] = 235.88
        v_['A3'] = 9.466
        v_['B3'] = -3347.25
        v_['C3'] = 215.31
        v_['AL1'] = 0.0
        v_['ALp1'] = 37.6
        v_['ALpp1'] = 0.0
        v_['AL2'] = 0.0
        v_['ALp2'] = 48.2
        v_['ALpp2'] = 0.0
        v_['AL3'] = 0.0
        v_['ALp3'] = 45.4
        v_['ALpp3'] = 0.0
        v_['BE1'] = 8425.0
        v_['BEp1'] = 24.2
        v_['BEpp1'] = 0.0
        v_['BE2'] = 9395.0
        v_['BEp2'] = 35.6
        v_['BEpp2'] = 0.0
        v_['BE3'] = 10466.0
        v_['BEp3'] = 31.9
        v_['BEpp3'] = 0.0
        v_['FL1'] = 30.0
        v_['FL2'] = 30.0
        v_['FL3'] = 40.0
        v_['FV1'] = 0.0
        v_['FV2'] = 0.0
        v_['FV3'] = 0.0
        v_['TF'] = 100.0
        v_['B'] = 40.0
        v_['D'] = 60.0
        v_['Q'] = 2500000.0
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['PI'+str(I)] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
            v_['INVPI'+str(I)] = 1.0/v_['PI'+str(I)]
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        for I in range(int(v_['0']),int(v_['N-2'])+1):
            [iv,ix_,_] = s2mpj_ii('V'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('2.1-'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'2.1-'+str(J))
            iv = ix_['X'+str(int(v_['0']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['B'])+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+2))
            [ig,ig_,_] = s2mpj_ii('2.3-'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'2.3-'+str(J))
            iv = ix_['X'+str(int(v_['N-1']))+','+str(J)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                [ig,ig_,_] = s2mpj_ii('2.2-'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'2.2-'+str(I)+','+str(J))
                pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+2))
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('2.7-'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'2.7-'+str(I))
        [ig,ig_,_] = s2mpj_ii('2.8',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'2.8')
        pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+5))
        for I in range(int(v_['1']),int(v_['N-2'])+1):
            [ig,ig_,_] = s2mpj_ii('2.9-'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'2.9-'+str(I))
            pbm.gscale = arrset(pbm.gscale,ig,float(1.0e+5))
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
        v_['SMALLHF'] = 0.0e+0
        v_['BIGHF'] = 0.0e+0
        for J in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['2.2-'+str(int(v_['K']))+','+str(J)],float(v_['FL'+str(J)])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['2.2-'+str(int(v_['K+']))+','+str(J)],float(v_['FV'+str(J)])))
            v_['TFTF'] = v_['TF']*v_['TF']
            v_['TEMP1'] = v_['TFTF']*v_['ALpp'+str(J)]
            v_['TEMP2'] = v_['TF']*v_['ALp'+str(J)]
            v_['TEMP1'] = v_['TEMP1']+v_['TEMP2']
            v_['TEMP1'] = v_['TEMP1']+v_['AL'+str(J)]
            v_['TEMP1'] = v_['TEMP1']*v_['FL'+str(J)]
            v_['SMALLHF'] = v_['SMALLHF']+v_['TEMP1']
            v_['TEMP1'] = v_['TFTF']*v_['BEpp'+str(J)]
            v_['TEMP2'] = v_['TF']*v_['BEp'+str(J)]
            v_['TEMP1'] = v_['TEMP1']+v_['TEMP2']
            v_['TEMP1'] = v_['TEMP1']+v_['BE'+str(J)]
            v_['TEMP1'] = v_['TEMP1']*v_['FV'+str(J)]
            v_['BIGHF'] = v_['BIGHF']+v_['TEMP1']
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['2.7-'+str(I)],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['2.8'],float(v_['Q']))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['2.9-'+str(int(v_['K']))],float(v_['SMALLHF'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['2.9-'+str(int(v_['K+']))],float(v_['BIGHF'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X0,1' in ix_):
            pb.x0[ix_['X0,1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X0,1']),float(0.0)))
        if('X0,2' in ix_):
            pb.x0[ix_['X0,2']] = float(0.2)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X0,2']),float(0.2)))
        if('X0,3' in ix_):
            pb.x0[ix_['X0,3']] = float(0.9)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X0,3']),float(0.9)))
        if('X1,1' in ix_):
            pb.x0[ix_['X1,1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1,1']),float(0.0)))
        if('X1,2' in ix_):
            pb.x0[ix_['X1,2']] = float(0.2)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1,2']),float(0.2)))
        if('X1,3' in ix_):
            pb.x0[ix_['X1,3']] = float(0.8)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1,3']),float(0.8)))
        if('X2,1' in ix_):
            pb.x0[ix_['X2,1']] = float(0.05)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2,1']),float(0.05)))
        if('X2,2' in ix_):
            pb.x0[ix_['X2,2']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2,2']),float(0.3)))
        if('X2,3' in ix_):
            pb.x0[ix_['X2,3']] = float(0.8)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2,3']),float(0.8)))
        if('X3,1' in ix_):
            pb.x0[ix_['X3,1']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3,1']),float(0.1)))
        if('X3,2' in ix_):
            pb.x0[ix_['X3,2']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3,2']),float(0.3)))
        if('X3,3' in ix_):
            pb.x0[ix_['X3,3']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3,3']),float(0.6)))
        if('X4,1' in ix_):
            pb.x0[ix_['X4,1']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,1']),float(0.3)))
        if('X4,2' in ix_):
            pb.x0[ix_['X4,2']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,2']),float(0.5)))
        if('X4,3' in ix_):
            pb.x0[ix_['X4,3']] = float(0.3)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4,3']),float(0.3)))
        if('X5,1' in ix_):
            pb.x0[ix_['X5,1']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,1']),float(0.6)))
        if('X5,2' in ix_):
            pb.x0[ix_['X5,2']] = float(0.6)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,2']),float(0.6)))
        if('X5,3' in ix_):
            pb.x0[ix_['X5,3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5,3']),float(0.0)))
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            if('T'+str(I) in ix_):
                pb.x0[ix_['T'+str(I)]] = float(100.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['T'+str(I)]),float(100.0)))
        for I in range(int(v_['0']),int(v_['N-2'])+1):
            if('V'+str(I) in ix_):
                pb.x0[ix_['V'+str(I)]] = float(300.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['V'+str(I)]),float(300.0)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        [it,iet_,_] = s2mpj_ii( 'ePOLY1PRD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P6')
        elftp = loaset(elftp,it,2,'P7')
        elftp = loaset(elftp,it,3,'P8')
        [it,iet_,_] = s2mpj_ii( 'ePOLY2PRD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P6')
        elftp = loaset(elftp,it,3,'P7')
        elftp = loaset(elftp,it,4,'P8')
        [it,iet_,_] = s2mpj_ii( 'eEXP2PROD', iet_)
        elftv = loaset(elftv,it,0,'V2')
        elftv = loaset(elftv,it,1,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        [it,iet_,_] = s2mpj_ii( 'eEXP3PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        [it,iet_,_] = s2mpj_ii( 'eEXP4PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'P1')
        elftp = loaset(elftp,it,1,'P2')
        elftp = loaset(elftp,it,2,'P3')
        elftp = loaset(elftp,it,3,'P4')
        elftp = loaset(elftp,it,4,'P5')
        elftp = loaset(elftp,it,5,'P6')
        elftp = loaset(elftp,it,6,'P7')
        elftp = loaset(elftp,it,7,'P8')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        v_['-D'] = -1.0*v_['D']
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E11-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
            ielftype = arrset(ielftype, ie, iet_["en2PROD"])
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'V'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            ename = 'E12-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP3PROD')
            ielftype = arrset(ielftype, ie, iet_["eEXP3PROD"])
            vname = 'V'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['0']))])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                v_['I-1'] = -1+I
                v_['I+1'] = 1+I
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype, ie, iet_["en2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
                ename = 'E22-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXP3PROD')
                ielftype = arrset(ielftype, ie, iet_["eEXP3PROD"])
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['I-1']))])))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype, ie, iet_["en2PROD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                ename = 'E24-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXP3PROD')
                ielftype = arrset(ielftype, ie, iet_["eEXP3PROD"])
                vname = 'V'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
            for I in range(int(v_['1']),int(v_['K-'])+1):
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            ename = 'E21-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
            ename = 'E23-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            for I in range(int(v_['K+']),int(v_['N-2'])+1):
                ename = 'E21-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
                ename = 'E23-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
            ename = 'E31-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP2PROD')
            ielftype = arrset(ielftype, ie, iet_["eEXP2PROD"])
            vname = 'X'+str(int(v_['N-2']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['N-2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['N-2']))])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['0']),int(v_['N-1'])+1):
                ename = 'E71-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXP2PROD')
                ielftype = arrset(ielftype, ie, iet_["eEXP2PROD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E81-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP4PROD')
            ielftype = arrset(ielftype, ie, iet_["eEXP4PROD"])
            vname = 'V'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['0']))])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BE'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
            ename = 'E82-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePOLY1PRD')
            ielftype = arrset(ielftype, ie, iet_["ePOLY1PRD"])
            vname = 'X'+str(int(v_['0']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['AL'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            ename = 'E83-'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePOLY2PRD')
            ielftype = arrset(ielftype, ie, iet_["ePOLY2PRD"])
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'V'+str(int(v_['0']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['AL'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                v_['I-1'] = -1+I
                v_['I+1'] = 1+I
                ename = 'E91-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXP4PROD')
                ielftype = arrset(ielftype, ie, iet_["eEXP4PROD"])
                vname = 'V'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(I)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BE'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePOLY2PRD')
                ielftype = arrset(ielftype, ie, iet_["ePOLY2PRD"])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['AL'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
                ename = 'E93-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXP4PROD')
                ielftype = arrset(ielftype, ie, iet_["eEXP4PROD"])
                vname = 'V'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I-1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['INVPI'+str(int(v_['I-1']))])))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P3')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P4')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P5')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BE'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEp'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BEpp'+str(J)]))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePOLY2PRD')
                ielftype = arrset(ielftype, ie, iet_["ePOLY2PRD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'T'+str(int(v_['I+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P6')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['AL'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P7')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALp'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='P8')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALpp'+str(J)]))
            for I in range(int(v_['1']),int(v_['K-'])+1):
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            ename = 'E92-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B']))
            ename = 'E94-'+str(int(v_['K']))+','+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
            for I in range(int(v_['K+']),int(v_['N-2'])+1):
                ename = 'E92-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
                ename = 'E94-'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='P2')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-D']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['2.1-'+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E11-'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E12-'+str(J)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_['2.3-'+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E31-'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['2.8']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E81-'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E82-'+str(J)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E83-'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            for I in range(int(v_['1']),int(v_['N-2'])+1):
                ig = ig_['2.2-'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E21-'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E22-'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E23-'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E24-'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                ig = ig_['2.9-'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E91-'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E92-'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E93-'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E94-'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            for I in range(int(v_['0']),int(v_['N-1'])+1):
                ig = ig_['2.7-'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E71-'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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
        pb.pbclass = "NOR2-AN-29-29"
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
        f_   = pbm.elpar[iel_][0]*EV_[0]*(EV_[1]+pbm.elpar[iel_][1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*(EV_[1]+pbm.elpar[iel_][1])
            g_[1] = pbm.elpar[iel_][0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOLY1PRD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        POLY  = (
              pbm.elpar[iel_][1]+pbm.elpar[iel_][2]*EV_[1]+pbm.elpar[iel_][3]*EV_[1]*EV_[1])
        DPOLY = pbm.elpar[iel_][2]+2.0*pbm.elpar[iel_][3]*EV_[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*POLY
            g_[1] = pbm.elpar[iel_][0]*EV_[0]*DPOLY
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]*DPOLY
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.elpar[iel_][0]*EV_[0]*2.0e+0*pbm.elpar[iel_][3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePOLY2PRD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        POLY  = (
              pbm.elpar[iel_][2]+pbm.elpar[iel_][3]*EV_[2]+pbm.elpar[iel_][4]*EV_[2]*EV_[2])
        DPOLY = pbm.elpar[iel_][3]+2.0*pbm.elpar[iel_][4]*EV_[2]
        f_   = pbm.elpar[iel_][0]*EV_[0]*(pbm.elpar[iel_][1]+EV_[1])*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*(pbm.elpar[iel_][1]+EV_[1])*POLY
            g_[1] = pbm.elpar[iel_][0]*EV_[0]*POLY
            g_[2] = pbm.elpar[iel_][0]*EV_[0]*(pbm.elpar[iel_][1]+EV_[1])*DPOLY
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = pbm.elpar[iel_][0]*POLY
                H_[1,0] = H_[0,1]
                H_[0,2] = pbm.elpar[iel_][0]*(pbm.elpar[iel_][1]+EV_[1])*DPOLY
                H_[2,0] = H_[0,2]
                H_[1,2] = pbm.elpar[iel_][0]*EV_[0]*DPOLY
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      pbm.elpar[iel_][0]*EV_[0]*(pbm.elpar[iel_][1]+EV_[1])*2.0e+0*pbm.elpar[iel_][4])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP2PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              pbm.elpar[iel_][0]*pbm.elpar[iel_][1]*np.exp(pbm.elpar[iel_][2]+(pbm.elpar[iel_][3]/(EV_[1]+pbm.elpar[iel_][4]))))
        F = EV_[0]*EXPROD
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPROD
            g_[1] = -EV_[0]*EXPROD*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[1])**2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -EXPROD*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[1])**2
                H_[1,0] = H_[0,1]
                H_[1,1] = (F*(pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[1])**2)**2+
                     2.0e+0*F*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[1])**3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP3PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              pbm.elpar[iel_][0]*pbm.elpar[iel_][1]*np.exp(pbm.elpar[iel_][2]+(pbm.elpar[iel_][3]/(EV_[2]+pbm.elpar[iel_][4]))))
        F = EV_[0]*EV_[1]*EXPROD
        TERM = -pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EXPROD
            g_[1] = EV_[0]*EXPROD
            g_[2] = F*TERM
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EXPROD
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EXPROD*TERM
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*EXPROD*TERM
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      F*(TERM*TERM+2.0e+0*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**3))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEXP4PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPROD  = (
              pbm.elpar[iel_][0]*pbm.elpar[iel_][1]*np.exp(pbm.elpar[iel_][2]+(pbm.elpar[iel_][3]/(EV_[2]+pbm.elpar[iel_][4]))))
        F = EV_[0]*EV_[1]*EXPROD
        POLY  = (
              pbm.elpar[iel_][5]+pbm.elpar[iel_][6]*EV_[2]+pbm.elpar[iel_][7]*EV_[2]*EV_[2])
        DPOLY = pbm.elpar[iel_][6]+2.0*pbm.elpar[iel_][7]*EV_[2]
        TERM = DPOLY-POLY*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**2
        f_   = F*POLY
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EXPROD*POLY
            g_[1] = EV_[0]*EXPROD*POLY
            g_[2] = F*TERM
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EXPROD*POLY
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EXPROD*TERM
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*EXPROD*TERM
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      F*(-(pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**2)*TERM+2.0*pbm.elpar[iel_][7]-DPOLY*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**2+2.0e+0*POLY*pbm.elpar[iel_][3]/(pbm.elpar[iel_][4]+EV_[2])**3))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

