from s2mpjlib import *
class  MRIBASIS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ******** 
# 
#    An optmization problem arising in the design of medical apparatus.
# 
#    Source:
#    Contribution from a LANCELOT user.
# 
#    SIF input: Arie Quist, TU Delft (NL), 1994.
#    Adaptation for CUTE: Ph. Toint, November 1994.
# 
#    classification = "LOR2-MY-36-55"
# 
#    useful constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MRIBASIS'

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
        v_['2'] = 2
        v_['3'] = 3
        v_['R2'] = float(v_['2'])
        v_['k1'] = 22.8443
        v_['k2'] = 12.4402
        v_['k3'] = 5.23792
        v_['k4'] = 5.12238
        v_['k5'] = 6.44999
        v_['k6'] = 5.32383
        v_['k7'] = 0.58392
        v_['k8'] = 3.94584
        v_['k9'] = -2.75584
        v_['k10'] = 32.0669
        v_['k11'] = 18.2179
        v_['k12'] = 31.7496
        v_['S1,1'] = -0.377126
        v_['S1,2'] = 0.919679
        v_['S1,3'] = 0.109389
        v_['S2,1'] = 0.634857
        v_['S2,2'] = 0.170704
        v_['S2,3'] = 0.753536
        v_['S3,1'] = 0.674338
        v_['S3,2'] = 0.353624
        v_['S3,3'] = -0.648242
        v_['xlo1'] = 1*v_['k7']
        v_['k10/2'] = v_['k10']/v_['R2']
        v_['k8/2'] = v_['k8']/v_['R2']
        v_['xup1'] = v_['k10/2']-v_['k8/2']
        v_['xlo2'] = v_['k10/2']+v_['k8/2']
        v_['k13'] = v_['k7']*v_['k5']
        v_['k14'] = v_['k8/2']*v_['k6']
        v_['xup2'] = 1*v_['k12']
        v_['-k1'] = -1*v_['k1']
        v_['-k2'] = -1*v_['k2']
        v_['k4-'] = v_['k4']-v_['k14']
        v_['-k3'] = -1*v_['k3']
        v_['-S1,1'] = -1*v_['S1,1']
        v_['-S1,2'] = -1*v_['S1,2']
        v_['-S1,3'] = -1*v_['S1,3']
        v_['-S2,1'] = -1*v_['S2,1']
        v_['-S2,2'] = -1*v_['S2,2']
        v_['-S2,3'] = -1*v_['S2,3']
        v_['-S3,1'] = -1*v_['S3,1']
        v_['-S3,2'] = -1*v_['S3,2']
        v_['-S3,3'] = -1*v_['S3,3']
        v_['2S1,1'] = 2*v_['S1,1']
        v_['2S1,2'] = 2*v_['S1,2']
        v_['2S1,3'] = 2*v_['S1,3']
        v_['2S2,1'] = 2*v_['S2,1']
        v_['2S2,2'] = 2*v_['S2,2']
        v_['2S2,3'] = 2*v_['S2,3']
        v_['2S3,1'] = 2*v_['S3,1']
        v_['2S3,2'] = 2*v_['S3,2']
        v_['2S3,3'] = 2*v_['S3,3']
        v_['-2S1,1'] = -2*v_['S1,1']
        v_['-2S1,2'] = -2*v_['S1,2']
        v_['-2S1,3'] = -2*v_['S1,3']
        v_['-2S2,1'] = -2*v_['S2,1']
        v_['-2S2,2'] = -2*v_['S2,2']
        v_['-2S2,3'] = -2*v_['S2,3']
        v_['-2S3,1'] = -2*v_['S3,1']
        v_['-2S3,2'] = -2*v_['S3,2']
        v_['-2S3,3'] = -2*v_['S3,3']
        v_['Llo1,1'] = v_['S1,1']*v_['k5']
        v_['Llo1,2'] = v_['S1,2']*v_['k5']
        v_['Llo1,3'] = v_['S1,3']*v_['k5']
        v_['Lup1,1'] = v_['S1,1']*v_['k6']
        v_['Lup1,2'] = v_['S1,2']*v_['k6']
        v_['Lup1,3'] = v_['S1,3']*v_['k6']
        v_['Llo2,1'] = v_['S1,1']*v_['k6']
        v_['Llo2,2'] = v_['S1,2']*v_['k6']
        v_['Llo2,3'] = v_['S1,3']*v_['k6']
        v_['4'] = 4
        v_['xm'] = 6
        v_['Lm'] = 4
        v_['xm-'] = -1+v_['xm']
        v_['xm-2'] = -2+v_['xm']
        v_['Lm-'] = -1+v_['Lm']
        v_['Lm-2'] = -2+v_['Lm']
        v_['R12'] = 12
        v_['1/12'] = 1/v_['R12']
        v_['-1/12'] = -1*v_['1/12']
        v_['R0'] = float(v_['0'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for j in range(int(v_['1']),int(v_['2'])+1):
            for k in range(int(v_['1']),int(v_['xm'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(j)+','+str(k),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(j)+','+str(k))
        for i in range(int(v_['1']),int(v_['3'])+1):
            for j in range(int(v_['1']),int(v_['2'])+1):
                for k in range(int(v_['1']),int(v_['Lm'])+1):
                    [iv,ix_,_] = s2mpj_ii('L'+str(i)+','+str(j)+','+str(k),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'L'+str(i)+','+str(j)+','+str(k))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('Object',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X'+str(int(v_['2']))+','+str(int(v_['xm']))]
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        for j in range(int(v_['1']),int(v_['2'])+1):
            for k in range(int(v_['2']),int(v_['xm-2'])+1):
                v_['k+'] = 1+k
                [ig,ig_,_] = s2mpj_ii('PS'+str(j)+','+str(k),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PS'+str(j)+','+str(k))
                iv = ix_['X'+str(j)+','+str(int(v_['k+']))]
                pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
                iv = ix_['X'+str(j)+','+str(k)]
                pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('PL',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'PL')
        iv = ix_['X'+str(int(v_['2']))+','+str(int(v_['xm']))]
        pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
        iv = ix_['X'+str(int(v_['2']))+','+str(int(v_['xm-']))]
        pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['3'])+1):
            for j in range(int(v_['1']),int(v_['2'])+1):
                for k in range(int(v_['1']),int(v_['Lm-'])+1):
                    v_['k+'] = 1+k
                    v_['2k'] = 2*k
                    v_['2k-'] = -1+v_['2k']
                    [ig,ig_,_] = s2mpj_ii('SU'+str(i)+','+str(j)+','+str(k),ig_)
                    gtype = arrset(gtype,ig,'<=')
                    cnames = arrset(cnames,ig,'SU'+str(i)+','+str(j)+','+str(k))
                    iv = ix_['L'+str(i)+','+str(j)+','+str(int(v_['k+']))]
                    pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
                    iv = ix_['L'+str(i)+','+str(j)+','+str(k)]
                    pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
                    iv = ix_['X'+str(j)+','+str(int(v_['2k']))]
                    pbm.A[ig,iv] = float(v_['-k1'])+pbm.A[ig,iv]
                    iv = ix_['X'+str(j)+','+str(int(v_['2k-']))]
                    pbm.A[ig,iv] = float(v_['k1'])+pbm.A[ig,iv]
                    [ig,ig_,_] = s2mpj_ii('SL'+str(i)+','+str(j)+','+str(k),ig_)
                    gtype = arrset(gtype,ig,'>=')
                    cnames = arrset(cnames,ig,'SL'+str(i)+','+str(j)+','+str(k))
                    iv = ix_['L'+str(i)+','+str(j)+','+str(int(v_['k+']))]
                    pbm.A[ig,iv] = float(1)+pbm.A[ig,iv]
                    iv = ix_['L'+str(i)+','+str(j)+','+str(k)]
                    pbm.A[ig,iv] = float(-1)+pbm.A[ig,iv]
                    iv = ix_['X'+str(j)+','+str(int(v_['2k']))]
                    pbm.A[ig,iv] = float(v_['k1'])+pbm.A[ig,iv]
                    iv = ix_['X'+str(j)+','+str(int(v_['2k-']))]
                    pbm.A[ig,iv] = float(v_['-k1'])+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['3'])+1):
            for k in range(int(v_['2']),int(v_['Lm-'])+1):
                [ig,ig_,_] = s2mpj_ii('cc1',ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'cc1')
                iv = ix_['L'+str(i)+','+str(int(v_['2']))+','+str(k)]
                pbm.A[ig,iv] = float(v_['S'+str(int(v_['3']))+','+str(i)])+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['3'])+1):
            [ig,ig_,_] = s2mpj_ii('c2const'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'c2const'+str(i))
            iv = ix_['L'+str(i)+','+str(int(v_['2']))+','+str(int(v_['Lm']))]
            pbm.A[ig,iv] = float(v_['k10'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('c3con1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'c3con1')
        [ig,ig_,_] = s2mpj_ii('c3con2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'c3con2')
        [ig,ig_,_] = s2mpj_ii('c4const',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'c4const')
        for i in range(int(v_['1']),int(v_['3'])+1):
            [ig,ig_,_] = s2mpj_ii('c5con'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'c5con'+str(i))
        for i in range(int(v_['1']),int(v_['2'])+1):
            [ig,ig_,_] = s2mpj_ii('c6cn'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'c6cn'+str(i))
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
        v_['Opmr1'] = v_['k3']*v_['S2,1']
        v_['Opmr2'] = v_['k3']*v_['S2,2']
        v_['Opmr3'] = v_['k3']*v_['S2,3']
        for i in range(int(v_['1']),int(v_['3'])+1):
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['c2const'+str(i)],float(v_['Opmr'+str(i)])))
        pbm.gconst = arrset(pbm.gconst,ig_['c3con1'],float(v_['k4-']))
        pbm.gconst = arrset(pbm.gconst,ig_['c3con2'],float(v_['k13']))
        pbm.gconst = arrset(pbm.gconst,ig_['c5con1'],float(v_['k13']))
        pbm.gconst = arrset(pbm.gconst,ig_['c5con2'],float(v_['-k3']))
        pbm.gconst = arrset(pbm.gconst,ig_['c5con3'],float(v_['k9']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for j in range(int(v_['1']),int(v_['2'])+1):
            for k in range(int(v_['2']),int(v_['xm-'])+1):
                pb.xlower[ix_['X'+str(j)+','+str(k)]] = v_['xlo'+str(j)]
                pb.xupper[ix_['X'+str(j)+','+str(k)]] = v_['xup'+str(j)]
            pb.xlower[ix_['X'+str(j)+','+str(int(v_['1']))]] = v_['xlo'+str(j)]
            pb.xupper[ix_['X'+str(j)+','+str(int(v_['1']))]] = v_['xlo'+str(j)]
        pb.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['xm']))]] = v_['xup1']
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['xm']))]] = v_['xup1']
        pb.xlower[ix_['X'+str(int(v_['2']))+','+str(int(v_['xm']))]] = v_['k11']
        pb.xupper[ix_['X'+str(int(v_['2']))+','+str(int(v_['xm']))]] = v_['k12']
        for i in range(int(v_['1']),int(v_['3'])+1):
            for j in range(int(v_['1']),int(v_['2'])+1):
                for k in range(int(v_['2']),int(v_['Lm-'])+1):
                    pb.xlower[ix_['L'+str(i)+','+str(j)+','+str(k)]] = v_['-k2']
                    pb.xupper[ix_['L'+str(i)+','+str(j)+','+str(k)]] = v_['k2']
                pb.xlower[ix_['L'+str(i)+','+str(j)+','+str(int(v_['1']))]] = (v_['Llo'+
                     str(j)+','+str(i)])
                pb.xupper[ix_['L'+str(i)+','+str(j)+','+str(int(v_['1']))]] = (v_['Llo'+
                     str(j)+','+str(i)])
            pb.xlower[ix_['L'+str(i)+','+str(int(v_['1']))+','+str(int(v_['Lm']))]]  = (
                  v_['Lup'+str(int(v_['1']))+','+str(i)])
            pb.xupper[ix_['L'+str(i)+','+str(int(v_['1']))+','+str(int(v_['Lm']))]]  = (
                  v_['Lup'+str(int(v_['1']))+','+str(i)])
            pb.xlower[ix_['L'+str(i)+','+str(int(v_['2']))+','+str(int(v_['Lm']))]]  = (
                  v_['-k2'])
            pb.xupper[ix_['L'+str(i)+','+str(int(v_['2']))+','+str(int(v_['Lm']))]]  = (
                  v_['k2'])
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        v_['intlen1'] = v_['xup1']-v_['xlo1']
        v_['Rxm-'] = float(v_['xm-'])
        v_['dx1'] = v_['intlen1']/v_['Rxm-']
        v_['intlen2'] = v_['k11']-v_['xlo2']
        v_['dx2'] = v_['intlen2']/v_['Rxm-']
        for k in range(int(v_['1']),int(v_['xm-'])+1):
            v_['Rk'] = float(k)
            v_['dist1'] = v_['dx1']*v_['Rk']
            v_['strtv1'] = v_['xlo1']+v_['dist1']
            v_['dist2'] = v_['dx2']*v_['Rk']
            v_['strtv2'] = v_['xlo2']+v_['dist2']
            v_['k+'] = 1+k
            pb.x0[ix_['X'+str(int(v_['1']))+','+str(int(v_['k+']))]]  = (
                  float(v_['strtv1']))
            pb.x0[ix_['X'+str(int(v_['2']))+','+str(int(v_['k+']))]]  = (
                  float(v_['strtv2']))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'euv1', iet_)
        elftv = loaset(elftv,it,0,'v1')
        elftv = loaset(elftv,it,1,'v2')
        elftv = loaset(elftv,it,2,'v3')
        elftv = loaset(elftv,it,3,'v4')
        [it,iet_,_] = s2mpj_ii( 'euv2', iet_)
        elftv = loaset(elftv,it,0,'v1')
        elftv = loaset(elftv,it,1,'v2')
        elftv = loaset(elftv,it,2,'v3')
        [it,iet_,_] = s2mpj_ii( 'euvw1', iet_)
        elftv = loaset(elftv,it,0,'v1')
        elftv = loaset(elftv,it,1,'v2')
        elftv = loaset(elftv,it,2,'v3')
        elftp = []
        elftp = loaset(elftp,it,0,'p1')
        [it,iet_,_] = s2mpj_ii( 'emo', iet_)
        elftv = loaset(elftv,it,0,'s1')
        elftv = loaset(elftv,it,1,'s2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for i in range(int(v_['1']),int(v_['3'])+1):
            for j in range(int(v_['1']),int(v_['2'])+1):
                for k in range(int(v_['1']),int(v_['Lm-'])+1):
                    v_['2k'] = 2*k
                    v_['k+'] = 1+k
                    v_['2k-'] = -1+v_['2k']
                    ename = 'e1'+str(i)+','+str(j)+','+str(k)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'euv1')
                    ielftype = arrset(ielftype, ie, iet_["euv1"])
                    vname = 'X'+str(j)+','+str(int(v_['2k']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k-']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'L'+str(i)+','+str(j)+','+str(int(v_['k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'L'+str(i)+','+str(j)+','+str(k)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v4')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'e3'+str(i)+','+str(j)+','+str(k)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'euvw1')
                    ielftype = arrset(ielftype, ie, iet_["euvw1"])
                    vname = 'L'+str(i)+','+str(j)+','+str(k)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k-']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    posep = find(elftp[ielftype[ie]],lambda x:x=='p1')
                    pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['1/12']))
                    ename = 'e5'+str(i)+','+str(j)+','+str(k)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'euvw1')
                    ielftype = arrset(ielftype, ie, iet_["euvw1"])
                    vname = 'L'+str(i)+','+str(j)+','+str(int(v_['k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k-']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    posep = find(elftp[ielftype[ie]],lambda x:x=='p1')
                    pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['-1/12']))
                for k in range(int(v_['1']),int(v_['Lm-2'])+1):
                    v_['2k'] = 2*k
                    v_['k+'] = 1+k
                    v_['2k+'] = 1+v_['2k']
                    ename = 'e2'+str(i)+','+str(j)+','+str(k)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'euv2')
                    ielftype = arrset(ielftype, ie, iet_["euv2"])
                    vname = 'X'+str(j)+','+str(int(v_['2k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'L'+str(i)+','+str(j)+','+str(int(v_['k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'e4'+str(i)+','+str(j)+','+str(k)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'euvw1')
                    ielftype = arrset(ielftype, ie, iet_["euvw1"])
                    vname = 'L'+str(i)+','+str(j)+','+str(int(v_['k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k+']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(j)+','+str(int(v_['2k']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='v3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    posep = find(elftp[ielftype[ie]],lambda x:x=='p1')
                    pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['R0']))
        for i in range(int(v_['1']),int(v_['3'])+1):
            ename = 'factr'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'emo')
            ielftype = arrset(ielftype, ie, iet_["emo"])
            vname = 'X'+str(int(v_['2']))+','+str(int(v_['xm']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='s1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'L'+str(i)+','+str(int(v_['2']))+','+str(int(v_['Lm']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='s2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['1']),int(v_['3'])+1):
            ig = ig_['c2const'+str(i)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['factr'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for j in range(int(v_['1']),int(v_['2'])+1):
            for i in range(int(v_['1']),int(v_['3'])+1):
                for k in range(int(v_['1']),int(v_['Lm-'])+1):
                    ig = ig_['c3con'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e1'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['S'+str(int(v_['1']))+','+str(i)])))
                for k in range(int(v_['1']),int(v_['Lm-2'])+1):
                    ig = ig_['c3con'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e2'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['S'+str(int(v_['1']))+','+str(i)])))
        for i in range(int(v_['1']),int(v_['3'])+1):
            for k in range(int(v_['1']),int(v_['Lm-'])+1):
                ig = ig_['c4const']
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['e1'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw  = (
                      loaset(pbm.grelw,ig,posel,float(v_['S'+str(int(v_['2']))+','+str(i)])))
            for k in range(int(v_['1']),int(v_['Lm-2'])+1):
                ig = ig_['c4const']
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['e2'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw  = (
                      loaset(pbm.grelw,ig,posel,float(v_['S'+str(int(v_['2']))+','+str(i)])))
        for j in range(int(v_['1']),int(v_['3'])+1):
            for i in range(int(v_['1']),int(v_['3'])+1):
                for k in range(int(v_['1']),int(v_['Lm-'])+1):
                    ig = ig_['c5con'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e1'+str(i)+','+str(int(v_['1']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e1'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['S'+str(j)+','+str(i)]))
                for k in range(int(v_['1']),int(v_['Lm-2'])+1):
                    ig = ig_['c5con'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e2'+str(i)+','+str(int(v_['1']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e2'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['S'+str(j)+','+str(i)]))
        for j in range(int(v_['1']),int(v_['2'])+1):
            for i in range(int(v_['1']),int(v_['3'])+1):
                for k in range(int(v_['1']),int(v_['Lm-'])+1):
                    ig = ig_['c6cn'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e3'+str(i)+','+str(int(v_['1']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e5'+str(i)+','+str(int(v_['1']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e3'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e5'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-S'+str(j)+','+str(i)]))
                for k in range(int(v_['1']),int(v_['Lm-2'])+1):
                    ig = ig_['c6cn'+str(j)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e4'+str(i)+','+str(int(v_['1']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['2S'+str(j)+','+str(i)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['e4'+str(i)+','+str(int(v_['2']))+','+str(k)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-2S'+str(j)+','+str(i)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               18.2179000000
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
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
        pb.pbclass = "LOR2-MY-36-55"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def euv1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = 0.5e0*IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.5e0*IV_[1]
            g_[1] = 0.5e0*IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 0.5e0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def euv2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,2] = U_[1,2]+1
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
                H_[0,1] = 1.0e0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def euvw1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        d14 = 0.25e0
        f_    = (
              EV_[0]*(EV_[1]-EV_[2])*((d14+pbm.elpar[iel_][0])*EV_[1]+(d14-pbm.elpar[iel_][0])*EV_[2]))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  (EV_[1]-EV_[2])*((d14+pbm.elpar[iel_][0])*EV_[1]+(d14-pbm.elpar[iel_][0])*EV_[2]))
            g_[1]  = (
                  EV_[0]*2.0e0*(EV_[1]*(d14+pbm.elpar[iel_][0])-pbm.elpar[iel_][0]*EV_[2]))
            g_[2]  = (
                  EV_[0]*2.0e0*(-EV_[2]*(d14-pbm.elpar[iel_][0])-pbm.elpar[iel_][0]*EV_[1]))
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 2.0e0*(EV_[1]*(d14+pbm.elpar[iel_][0])-pbm.elpar[iel_][0]*EV_[2])
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0e0*(-EV_[2]*(d14-pbm.elpar[iel_][0])-pbm.elpar[iel_][0]*EV_[1])
                H_[2,0] = H_[0,2]
                H_[1,1] = EV_[0]*2.0e0*(d14+pbm.elpar[iel_][0])
                H_[1,2] = -EV_[0]*2.0e0*pbm.elpar[iel_][0]
                H_[2,1] = H_[1,2]
                H_[2,2] = -EV_[0]*2.0e0*(d14-pbm.elpar[iel_][0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def emo(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[1]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EV_[1]
            g_[1] = -EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -1.0e0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

