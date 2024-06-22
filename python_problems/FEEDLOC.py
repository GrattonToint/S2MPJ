from s2mpjlib import *
class  FEEDLOC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FEEDLOC
#    *********
# 
#    Feed tray location & determination of optimum number of trays 
#    in a distillation column
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "LOR2-AN-90-259"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FEEDLOC'

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
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['M'] = 2
        v_['NMAX'] = 12
        v_['NMAX-1'] = -1+v_['NMAX']
        v_['F'] = 100.0
        v_['AL1'] = 1.0
        v_['AL2'] = 5.13435
        v_['XF1'] = 0.80
        v_['XF2'] = 0.20
        v_['SPEC'] = 0.001
        v_['BIGM'] = 1000.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [iv,ix_,_] = s2mpj_ii('S'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'S'+str(I))
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I))
            [iv,ix_,_] = s2mpj_ii('Z'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Z'+str(I))
        [iv,ix_,_] = s2mpj_ii('N',ix_)
        pb.xnames=arrset(pb.xnames,iv,'N')
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(J))
        [iv,ix_,_] = s2mpj_ii('L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L')
        [iv,ix_,_] = s2mpj_ii('V',ix_)
        pb.xnames=arrset(pb.xnames,iv,'V')
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        pb.xnames=arrset(pb.xnames,iv,'R')
        [iv,ix_,_] = s2mpj_ii('P1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P1')
        [iv,ix_,_] = s2mpj_ii('P2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['R']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('FENTR',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'FENTR')
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('NTRAY',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NTRAY')
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('NDEF1',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NDEF1')
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            v_['RI'] = float(I)
            [ig,ig_,_] = s2mpj_ii('NDEF2',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NDEF2')
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(v_['RI'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('NDEF1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NDEF1')
        iv = ix_['N']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('NDEF2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'NDEF2')
        iv = ix_['N']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('NIL'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'NIL'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['Z'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            v_['RI'] = float(I)
            [ig,ig_,_] = s2mpj_ii('ENTX',ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'ENTX')
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['RI'])+pbm.A[ig,iv]
            v_['RI'] = -1.0*v_['RI']
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(v_['RI'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('LASTX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'LASTX'+str(I))
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('ZNOT'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'ZNOT'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            for K in range(int(I),int(v_['NMAX'])+1):
                [ig,ig_,_] = s2mpj_ii('ZNOT'+str(I),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'ZNOT'+str(I))
                iv = ix_['S'+str(K)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('FEEDX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'FEEDX'+str(I))
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('WNES1u'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'WNES1u'+str(I))
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('WNES2u'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'WNES2u'+str(I))
            iv = ix_['W'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('PE1'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE1'+str(I))
                iv = ix_['Y'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('PE2'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE2'+str(I))
                iv = ix_['Y'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('PE3'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE3'+str(I))
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('PE4'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'PE4'+str(I))
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('PE1'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE1'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('PE2'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE2'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('PE3'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE3'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('PE4'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'PE4'+str(I))
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('XNOT'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'XNOT'+str(I)+','+str(J))
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Z'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('YNOT'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'YNOT'+str(I)+','+str(J))
                iv = ix_['Y'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Z'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            v_['TEMP'] = -1.0*v_['AL1']
            [ig,ig_,_] = s2mpj_ii('PHEE'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PHEE'+str(I))
            iv = ix_['X'+str(I)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(v_['TEMP'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('DEFL',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'DEFL')
        iv = ix_['L']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['M'])+1):
            v_['TEMP'] = -1.0*v_['F']
            [ig,ig_,_] = s2mpj_ii('CMB1u'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CMB1u'+str(J))
            iv = ix_['X'+str(int(v_['2']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['TEMP'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                v_['TEMP'] = -1.0*v_['BIGM']
                [ig,ig_,_] = s2mpj_ii('CMBN1'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'CMBN1'+str(I)+','+str(J))
                iv = ix_['S'+str(I)]
                pbm.A[ig,iv] = float(v_['BIGM'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('CMBN2'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CMBN2'+str(I)+','+str(J))
                iv = ix_['S'+str(I)]
                pbm.A[ig,iv] = float(v_['TEMP'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            v_['TEMP1'] = v_['F']*v_['XF'+str(int(v_['M']))]
            v_['TEMP1'] = -1.0*v_['TEMP1']
            [ig,ig_,_] = s2mpj_ii('CMB1'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CMB1'+str(I))
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(v_['TEMP'])+pbm.A[ig,iv]
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(v_['BIGM'])+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['TEMP1'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('CMB2'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CMB2'+str(I))
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(v_['BIGM'])+pbm.A[ig,iv]
            iv = ix_['Z'+str(I)]
            pbm.A[ig,iv] = float(v_['TEMP'])+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['TEMP1'])+pbm.A[ig,iv]
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            [ig,ig_,_] = s2mpj_ii('RECR'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'RECR'+str(I))
            iv = ix_['S'+str(I)]
            pbm.A[ig,iv] = float(v_['BIGM'])+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['FENTR'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['NTRAY'],float(1.0))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['WNES1u'+str(I)],float(1.0))
            pbm.gconst = arrset(pbm.gconst,ig_['WNES2u'+str(I)],float(1.0))
        v_['TEMP'] = -1.0*v_['BIGM']
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['CMBN1'+str(I)+','+str(J)],float(v_['BIGM'])))
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['CMBN2'+str(I)+','+str(J)],float(v_['TEMP'])))
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CMB1'+str(I)],float(v_['BIGM']))
            v_['TEMP'] = -1.0*v_['BIGM']
            pbm.gconst = arrset(pbm.gconst,ig_['CMB2'+str(I)],float(v_['TEMP']))
        v_['TEMP'] = v_['XF'+str(int(v_['1']))]*v_['SPEC']
        v_['TEMP1'] = v_['TEMP']*v_['F']
        v_['RHS'] = v_['TEMP1']+v_['BIGM']
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['RECR'+str(I)],float(v_['RHS']))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            grange = arrset(grange,ig_['PE1'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE2'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE3'+str(I)],float(2.0))
            grange = arrset(grange,ig_['PE4'+str(I)],float(2.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            pb.xupper[ix_['Z'+str(I)]] = 1.0
            pb.xupper[ix_['W'+str(I)]] = 1.0
            pb.xupper[ix_['S'+str(I)]] = 1.0
            for J in range(int(v_['1']),int(v_['M'])+1):
                pb.xupper[ix_['X'+str(I)+','+str(J)]] = 1.0
                pb.xupper[ix_['Y'+str(I)+','+str(J)]] = 1.0
        pb.xlower[ix_['N']] = 3.0
        v_['TEMP'] = float(v_['NMAX'])
        pb.xupper[ix_['N']] = v_['TEMP']
        pb.xlower[ix_['P2']] = 80.0
        pb.xupper[ix_['P2']] = 80.0
        pb.xupper[ix_['L']] = v_['F']
        pb.xupper[ix_['V']] = v_['F']
        pb.xupper[ix_['P1']] = v_['F']
        pb.xupper[ix_['R']] = 5.0
        pb.xlower[ix_['W'+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['W'+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['W'+str(int(v_['2']))]] = 0.0
        pb.xupper[ix_['W'+str(int(v_['2']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        pb.y0 = np.full((pb.m,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eA2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                ename = 'PHE'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['AL'+str(K)]))
        ename = 'DEFLE'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
        ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
        vname = 'R'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'P1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='A')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'CMB11u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'P2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            ename = 'CMB12u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['1']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            ename = 'CMB13u'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['2']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['M'])+1):
                ename = 'CM11'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'L'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                ename = 'CM12'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'P1'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                ename = 'CM13'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'V'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
                ename = 'CM21'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'L'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                ename = 'CM22'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'P1'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
                ename = 'CM23'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'V'
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ename = 'C11'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            for K in range(int(I),int(v_['NMAX'])+1):
                ename = 'C12'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['F']))
            ename = 'C13'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            ename = 'C14'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                v_['TEMP'] = -1.0*v_['F']
                ename = 'C15'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TEMP']))
            ename = 'C16'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
            ename = 'C21'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            for K in range(int(v_['1']),int(v_['NMAX'])+1):
                ename = 'C22'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'X'+str(I)+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['F']))
            ename = 'C23'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            ename = 'C24'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'L'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                v_['TEMP'] = -1.0*v_['F']
                ename = 'C25'+str(I)+','+str(K)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
                ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
                vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'W'+str(K)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='A')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TEMP']))
            ename = 'C26'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'V'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I-1']))+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            ename = 'REC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eA2PROD')
            ielftype = arrset(ielftype, ie, iet_["eA2PROD"])
            vname = 'P1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)+','+str(int(v_['1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NMAX'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['PHEE'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PHE'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['DEFL']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DEFLE'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['CMB1u'+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CMB11u'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CMB12u'+str(J)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CMB13u'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['NMAX'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['CMBN1'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM11'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM12'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM13'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['CMBN2'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM21'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM22'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CM23'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['NMAX-1'])+1):
            ig = ig_['CMB1'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C11'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C13'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C14'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C16'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            ig = ig_['CMB2'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C21'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C23'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C24'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C26'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            for K in range(int(I),int(v_['NMAX'])+1):
                ig = ig_['CMB1'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C12'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['CMB2'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C22'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            v_['I+1'] = 1+I
            for K in range(int(v_['I+1']),int(v_['NMAX'])+1):
                ig = ig_['CMB1'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C15'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['CMB2'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C25'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['3']),int(v_['NMAX'])+1):
            ig = ig_['RECR'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['REC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
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
        pb.pbclass = "LOR2-AN-90-259"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eA2PROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EV_[1]
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

