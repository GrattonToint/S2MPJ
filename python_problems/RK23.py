from s2mpjlib import *
class  RK23(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RK23
#    *********
# 
#    Find coefficients for an embedded pair of explicit 2nd
#    and 3rd order Runge Kutta Method such that the leading
#    term in the local truncation error is minimized.
# 
#    Source:
#    Similar ideas for 4th and 5th order pairs are discussed in:
#    Hairer, Norsett and Wanner, Solving Ordinary Differential
#    Equations I, Springer 1980, page 158 ff.
# 
#    SIF input: S. Leyffer, January 1997.
# 
#    classification = "LOR2-RN-17-11"
# 
# 
#    ... COMPUTED PARAMETERS
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RK23'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'RK23'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['ONE'] = 1.0
        v_['THREE'] = 3.0
        v_['FOUR'] = 4.0
        v_['SIX'] = 6.0
        v_['ONETHIRD'] = v_['ONE']/v_['THREE']
        v_['ONESIXTH'] = v_['ONE']/v_['SIX']
        v_['FOURSIXTH'] = v_['FOUR']/v_['SIX']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('C2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C2')
        [iv,ix_,_] = s2mpj_ii('A21',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A21')
        [iv,ix_,_] = s2mpj_ii('C3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C3')
        [iv,ix_,_] = s2mpj_ii('A31',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A31')
        [iv,ix_,_] = s2mpj_ii('A32',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A32')
        [iv,ix_,_] = s2mpj_ii('B1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B1')
        [iv,ix_,_] = s2mpj_ii('B2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B2')
        [iv,ix_,_] = s2mpj_ii('B3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B3')
        [iv,ix_,_] = s2mpj_ii('BB1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BB1')
        [iv,ix_,_] = s2mpj_ii('BB2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BB2')
        [iv,ix_,_] = s2mpj_ii('BB3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'BB3')
        [iv,ix_,_] = s2mpj_ii('TP1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TP1')
        [iv,ix_,_] = s2mpj_ii('TM1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TM1')
        [iv,ix_,_] = s2mpj_ii('TP2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TP2')
        [iv,ix_,_] = s2mpj_ii('TM2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TM2')
        [iv,ix_,_] = s2mpj_ii('TP3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TP3')
        [iv,ix_,_] = s2mpj_ii('TM3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TM3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['TP1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TP2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TP3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('ROWS1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ROWS1')
        iv = ix_['A21']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['C2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('ROWS2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ROWS2')
        iv = ix_['A31']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['A32']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['C3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FIRST2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'FIRST2')
        iv = ix_['B1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['B2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['B3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('FIRST3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'FIRST3')
        iv = ix_['BB1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['BB2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['BB3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('SECND2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SECND2')
        [ig,ig_,_] = s2mpj_ii('SECND3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'SECND3')
        [ig,ig_,_] = s2mpj_ii('THIRD31',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'THIRD31')
        [ig,ig_,_] = s2mpj_ii('THIRD32',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'THIRD32')
        [ig,ig_,_] = s2mpj_ii('ART1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART1')
        iv = ix_['TP1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('ART2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART2')
        iv = ix_['TP2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('ART3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'ART3')
        iv = ix_['TP3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['TM3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['ROWS1'],float(0.0))
        pbm.gconst = arrset(pbm.gconst,ig_['ROWS2'],float(0.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FIRST2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['FIRST3'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['SECND2'],float(0.5))
        pbm.gconst = arrset(pbm.gconst,ig_['SECND3'],float(0.5))
        pbm.gconst = arrset(pbm.gconst,ig_['THIRD31'],float(v_['ONETHIRD']))
        pbm.gconst = arrset(pbm.gconst,ig_['THIRD32'],float(v_['ONESIXTH']))
        pbm.gconst = arrset(pbm.gconst,ig_['ART1'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['ART2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['ART3'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['TP1']] = 0.0
        pb.xlower[ix_['TM1']] = 0.0
        pb.xlower[ix_['TP2']] = 0.0
        pb.xlower[ix_['TM2']] = 0.0
        pb.xlower[ix_['TP3']] = 0.0
        pb.xlower[ix_['TM3']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('C2' in ix_):
            pb.x0[ix_['C2']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['C2']),float(1.0)))
        if('A21' in ix_):
            pb.x0[ix_['A21']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['A21']),float(1.0)))
        if('C3' in ix_):
            pb.x0[ix_['C3']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['C3']),float(0.5)))
        if('A31' in ix_):
            pb.x0[ix_['A31']] = float(0.25)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['A31']),float(0.25)))
        if('A32' in ix_):
            pb.x0[ix_['A32']] = float(0.25)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['A32']),float(0.25)))
        if('B1' in ix_):
            pb.x0[ix_['B1']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B1']),float(0.5)))
        if('B2' in ix_):
            pb.x0[ix_['B2']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B2']),float(0.5)))
        if('B3' in ix_):
            pb.x0[ix_['B3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B3']),float(0.0)))
        pb.x0[ix_['BB1']] = float(v_['ONESIXTH'])
        pb.x0[ix_['BB2']] = float(v_['ONESIXTH'])
        pb.x0[ix_['BB3']] = float(v_['FOURSIXTH'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'ePRODS', iet_)
        elftv = loaset(elftv,it,0,'W1')
        elftv = loaset(elftv,it,1,'W2')
        [it,iet_,_] = s2mpj_ii( 'ePRODQ', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        [it,iet_,_] = s2mpj_ii( 'eTPROD', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        [it,iet_,_] = s2mpj_ii( 'eQPROD', iet_)
        elftv = loaset(elftv,it,0,'Z1')
        elftv = loaset(elftv,it,1,'Z2')
        elftv = loaset(elftv,it,2,'Z3')
        elftv = loaset(elftv,it,3,'Z4')
        [it,iet_,_] = s2mpj_ii( 'eTPRODS', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'B2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'B3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'BB2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODS')
        ielftype = arrset(ielftype, ie, iet_["ePRODS"])
        vname = 'BB2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODS')
        ielftype = arrset(ielftype, ie, iet_["ePRODS"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='W2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eTPROD')
        ielftype = arrset(ielftype, ie, iet_["eTPROD"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODQ')
        ielftype = arrset(ielftype, ie, iet_["ePRODQ"])
        vname = 'BB2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E9'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePRODQ')
        ielftype = arrset(ielftype, ie, iet_["ePRODQ"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eQPROD')
        ielftype = arrset(ielftype, ie, iet_["eQPROD"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'E11'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eTPRODS')
        ielftype = arrset(ielftype, ie, iet_["eTPRODS"])
        vname = 'BB3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'A32'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'C2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['SECND2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['SECND3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['THIRD31']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['THIRD32']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['ART1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E9'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(4.0))
        ig = ig_['ART2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(8.0))
        ig = ig_['ART3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E11'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(12.0))
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
        pb.pbclass = "LOR2-RN-17-11"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*(EV_[1]**2.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]**2
            g_[1] = EV_[0]*2.0*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 2.0*EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePRODQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*(EV_[1]**3.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]**3.0
            g_[1] = EV_[0]*3.0*(EV_[1]**2.0)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 3.0*(EV_[1]**2.0)
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0]*6.0*EV_[1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTPROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQPROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]
            g_[1] = EV_[0]*EV_[2]*EV_[3]
            g_[2] = EV_[0]*EV_[1]*EV_[3]
            g_[3] = EV_[0]*EV_[1]*EV_[2]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2]*EV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0]*EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[0]*EV_[1]
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eTPRODS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*(EV_[2]**2.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*(EV_[2]**2.0)
            g_[1] = EV_[0]*(EV_[2]**2.0)
            g_[2] = EV_[0]*EV_[1]*2.0*EV_[2]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]**2.0
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*2.0*EV_[2]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]*2.0*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[0]*EV_[1]*2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

