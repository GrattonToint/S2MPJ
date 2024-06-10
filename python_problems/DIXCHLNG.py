from s2mpjlib import *
class  DIXCHLNG(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DIXCHLNG
#    *********
# 
#    A constrained problem set as a challenge for SQP methods
#    by L.C.W. Dixon at the APMOD91 Conference.
# 
#    Source:
#    L.C.W. Dixon, personnal communication, Jan 1991.
# 
#    SIF input: Ph. Toint, Feb 1991.
# 
#    classification = "SOR2-AN-10-5"
# 
#    Other parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DIXCHLNG'

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
        v_['7'] = 7
        v_['9'] = 9
        v_['10'] = 10
        v_['90.0'] = 90.0
        v_['10.1'] = 10.1
        v_['19.8'] = 19.8
        v_['1/90.0'] = 1.0/v_['90.0']
        v_['1/10.1'] = 1.0/v_['10.1']
        v_['1/19.8'] = 1.0/v_['19.8']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['10'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['I+3'] = 3+I
            [ig,ig_,_] = s2mpj_ii('A'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(0.01))
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+3']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/90.0']))
            [ig,ig_,_] = s2mpj_ii('D'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+2']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/10.1']))
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['I+3']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/10.1']))
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/19.8']))
        for I in range(int(v_['2']),int(v_['10'])+1,int(v_['2'])):
            [ig,ig_,_] = s2mpj_ii('P'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'P'+str(I))
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),1.0)
        for I in range(int(v_['1']),int(v_['7'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['A'+str(I)],float(0.0))
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)],float(0.0))
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        v_['X0A'] = 2.0
        v_['X0M'] = -1.0
        for I in range(int(v_['1']),int(v_['9'])+1,int(v_['2'])):
            v_['X0'] = v_['X0A']*v_['X0M']
            pb.x0[ix_['X'+str(I)]] = float(v_['X0'])
            v_['1/X0'] = 1.0/v_['X0']
            v_['I+1'] = 1+I
            pb.x0[ix_['X'+str(int(v_['I+1']))]] = float(v_['1/X0'])
            v_['X0A'] = 1.0+v_['X0A']
            v_['X0M'] = -1.0*v_['X0M']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eS2PR', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        [it,iet_,_] = s2mpj_ii( 'ePR2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'ePR4', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        [it,iet_,_] = s2mpj_ii( 'ePR6', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        [it,iet_,_] = s2mpj_ii( 'ePR8', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        [it,iet_,_] = s2mpj_ii( 'ePR10', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['9'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+1'] = 1+I
            v_['I+3'] = 3+I
            ename = 'PR'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eS2PR')
            ielftype = arrset(ielftype, ie, iet_["eS2PR"])
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PRD2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePR2')
        ielftype = arrset(ielftype, ie, iet_["ePR2"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PRD4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePR4')
        ielftype = arrset(ielftype, ie, iet_["ePR4"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PRD6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePR6')
        ielftype = arrset(ielftype, ie, iet_["ePR6"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PRD8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePR8')
        ielftype = arrset(ielftype, ie, iet_["ePR8"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V8')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'PRD10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePR10')
        ielftype = arrset(ielftype, ie, iet_["ePR10"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V8')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V9')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V10')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+2'] = 2+I
            ig = ig_['A'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['B'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['C'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(int(v_['I+2']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['D'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['F'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PR'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['10'])+1,int(v_['2'])):
            ig = ig_['P'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PRD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "SOR2-AN-10-5"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eS2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-1.0)*(EV_[1]-1.0)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]-1.0
            g_[1] = EV_[0]-1.0
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
    def ePR2(pbm,nargout,*args):

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
    def ePR4(pbm,nargout,*args):

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
    def ePR6(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]
            g_[1] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]
            g_[2] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]
            g_[3] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]
            g_[4] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]
            g_[5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[0,1] = EV_[2]*EV_[3]*EV_[4]*EV_[5]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]*EV_[4]*EV_[5]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]*EV_[4]*EV_[5]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1]*EV_[2]*EV_[3]*EV_[5]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1]*EV_[2]*EV_[3]*EV_[4]
                H_[5,0] = H_[0,5]
                H_[1,2] = EV_[0]*EV_[3]*EV_[4]*EV_[5]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]*EV_[4]*EV_[5]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0]*EV_[2]*EV_[3]*EV_[5]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0]*EV_[2]*EV_[3]*EV_[4]
                H_[5,1] = H_[1,5]
                H_[2,3] = EV_[0]*EV_[1]*EV_[4]*EV_[5]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0]*EV_[1]*EV_[3]*EV_[5]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0]*EV_[1]*EV_[3]*EV_[4]
                H_[5,2] = H_[2,5]
                H_[3,4] = EV_[0]*EV_[1]*EV_[2]*EV_[5]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0]*EV_[1]*EV_[2]*EV_[4]
                H_[5,3] = H_[3,5]
                H_[4,5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]
                H_[5,4] = H_[4,5]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePR8(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
            g_[1] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
            g_[2] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
            g_[3] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
            g_[4] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]
            g_[5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]
            g_[6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]
            g_[7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
            if nargout>2:
                H_ = np.zeros((8,8))
                H_[0,1] = EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]
                H_[5,0] = H_[0,5]
                H_[0,6] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]
                H_[6,0] = H_[0,6]
                H_[0,7] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
                H_[7,0] = H_[0,7]
                H_[1,2] = EV_[0]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]
                H_[6,1] = H_[1,6]
                H_[1,7] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
                H_[7,1] = H_[1,7]
                H_[2,3] = EV_[0]*EV_[1]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0]*EV_[1]*EV_[3]*EV_[5]*EV_[6]*EV_[7]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[6]*EV_[7]
                H_[5,2] = H_[2,5]
                H_[2,6] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[7]
                H_[6,2] = H_[2,6]
                H_[2,7] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]
                H_[7,2] = H_[2,7]
                H_[3,4] = EV_[0]*EV_[1]*EV_[2]*EV_[5]*EV_[6]*EV_[7]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[6]*EV_[7]
                H_[5,3] = H_[3,5]
                H_[3,6] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[7]
                H_[6,3] = H_[3,6]
                H_[3,7] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]
                H_[7,3] = H_[3,7]
                H_[4,5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[6]*EV_[7]
                H_[5,4] = H_[4,5]
                H_[4,6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[7]
                H_[6,4] = H_[4,6]
                H_[4,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]
                H_[7,4] = H_[4,7]
                H_[5,6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[7]
                H_[6,5] = H_[5,6]
                H_[5,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]
                H_[7,5] = H_[5,7]
                H_[6,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]
                H_[7,6] = H_[6,7]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePR10(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[1] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[2] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[3] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[4] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
            g_[6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
            g_[7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
            g_[8] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
            g_[9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,1] = EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[5,0] = H_[0,5]
                H_[0,6] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
                H_[6,0] = H_[0,6]
                H_[0,7] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
                H_[7,0] = H_[0,7]
                H_[0,8] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
                H_[8,0] = H_[0,8]
                H_[0,9] = EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
                H_[9,0] = H_[0,9]
                H_[1,2] = EV_[0]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
                H_[6,1] = H_[1,6]
                H_[1,7] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
                H_[7,1] = H_[1,7]
                H_[1,8] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
                H_[8,1] = H_[1,8]
                H_[1,9] = EV_[0]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
                H_[9,1] = H_[1,9]
                H_[2,3] = EV_[0]*EV_[1]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0]*EV_[1]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[5,2] = H_[2,5]
                H_[2,6] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
                H_[6,2] = H_[2,6]
                H_[2,7] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
                H_[7,2] = H_[2,7]
                H_[2,8] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
                H_[8,2] = H_[2,8]
                H_[2,9] = EV_[0]*EV_[1]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
                H_[9,2] = H_[2,9]
                H_[3,4] = EV_[0]*EV_[1]*EV_[2]*EV_[5]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[5,3] = H_[3,5]
                H_[3,6] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
                H_[6,3] = H_[3,6]
                H_[3,7] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
                H_[7,3] = H_[3,7]
                H_[3,8] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
                H_[8,3] = H_[3,8]
                H_[3,9] = EV_[0]*EV_[1]*EV_[2]*EV_[4]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
                H_[9,3] = H_[3,9]
                H_[4,5] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[6]*EV_[7]*EV_[8]*EV_[9]
                H_[5,4] = H_[4,5]
                H_[4,6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[7]*EV_[8]*EV_[9]
                H_[6,4] = H_[4,6]
                H_[4,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[8]*EV_[9]
                H_[7,4] = H_[4,7]
                H_[4,8] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[9]
                H_[8,4] = H_[4,8]
                H_[4,9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[5]*EV_[6]*EV_[7]*EV_[8]
                H_[9,4] = H_[4,9]
                H_[5,6] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[7]*EV_[8]*EV_[9]
                H_[6,5] = H_[5,6]
                H_[5,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[8]*EV_[9]
                H_[7,5] = H_[5,7]
                H_[5,8] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[9]
                H_[8,5] = H_[5,8]
                H_[5,9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[6]*EV_[7]*EV_[8]
                H_[9,5] = H_[5,9]
                H_[6,7] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[8]*EV_[9]
                H_[7,6] = H_[6,7]
                H_[6,8] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[9]
                H_[8,6] = H_[6,8]
                H_[6,9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[7]*EV_[8]
                H_[9,6] = H_[6,9]
                H_[7,8] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[9]
                H_[8,7] = H_[7,8]
                H_[7,9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[8]
                H_[9,7] = H_[7,9]
                H_[8,9] = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]*EV_[5]*EV_[6]*EV_[7]
                H_[9,8] = H_[8,9]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

