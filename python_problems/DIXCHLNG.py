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
#    classification = "C-CSOR2-AN-10-5"
# 
#    Other parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DIXCHLNG'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['10'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['I+3'] = 3+I
            [ig,ig_,_] = s2mpj_ii('A'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            self.gscale = arrset(self.gscale,ig,float(0.01))
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+3']))]])
            valA = np.append(valA,float(1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['1/90.0']))
            [ig,ig_,_] = s2mpj_ii('D'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+2']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['1/10.1']))
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+3']))]])
            valA = np.append(valA,float(1.0))
            self.gscale = arrset(self.gscale,ig,float(v_['1/10.1']))
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(v_['1/19.8']))
        for I in range(int(v_['2']),int(v_['10'])+1,int(v_['2'])):
            [ig,ig_,_] = s2mpj_ii('P'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'P'+str(I))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),1.0)
        for I in range(int(v_['1']),int(v_['7'])+1):
            self.gconst = arrset(self.gconst,ig_['A'+str(I)],float(0.0))
            self.gconst = arrset(self.gconst,ig_['C'+str(I)],float(0.0))
            self.gconst = arrset(self.gconst,ig_['G'+str(I)],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        v_['X0A'] = 2.0
        v_['X0M'] = -1.0
        for I in range(int(v_['1']),int(v_['9'])+1,int(v_['2'])):
            v_['X0'] = v_['X0A']*v_['X0M']
            self.x0[ix_['X'+str(I)]] = float(v_['X0'])
            v_['1/X0'] = 1.0/v_['X0']
            v_['I+1'] = 1+I
            self.x0[ix_['X'+str(int(v_['I+1']))]] = float(v_['1/X0'])
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['9'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+1'] = 1+I
            v_['I+3'] = 3+I
            ename = 'PR'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eS2PR')
            ielftype = arrset(ielftype,ie,iet_["eS2PR"])
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='W')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PRD2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePR2')
        ielftype = arrset(ielftype,ie,iet_["ePR2"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PRD4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePR4')
        ielftype = arrset(ielftype,ie,iet_["ePR4"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PRD6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePR6')
        ielftype = arrset(ielftype,ie,iet_["ePR6"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PRD8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePR8')
        ielftype = arrset(ielftype,ie,iet_["ePR8"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'PRD10'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePR10')
        ielftype = arrset(ielftype,ie,iet_["ePR10"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V4')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V5')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V6')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V7')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V8')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V9')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V10')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+2'] = 2+I
            ig = ig_['A'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['B'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['C'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['XSQ'+str(int(v_['I+2']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['D'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['F'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gL2')
            ig = ig_['G'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['PR'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['10'])+1,int(v_['2'])):
            ig = ig_['P'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['PRD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSOR2-AN-10-5"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[0,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0,0]+EV_[0,0]
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
    def eS2PR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0,0]-1.0)*(EV_[1,0]-1.0)
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]-1.0
            g_[1] = EV_[0,0]-1.0
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
    def ePR2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
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
    def ePR4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*EV_[3,0]
            g_[1] = EV_[0,0]*EV_[2,0]*EV_[3,0]
            g_[2] = EV_[0,0]*EV_[1,0]*EV_[3,0]
            g_[3] = EV_[0,0]*EV_[1,0]*EV_[2,0]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2,0]*EV_[3,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*EV_[3,0]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1,0]*EV_[2,0]
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0,0]*EV_[3,0]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0,0]*EV_[2,0]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[0,0]*EV_[1,0]
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePR6(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[1] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[2] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
            g_[3] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
            g_[4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
            g_[5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
            if nargout>2:
                H_ = np.zeros((6,6))
                H_[0,1] = EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
                H_[5,0] = H_[0,5]
                H_[1,2] = EV_[0,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]
                H_[5,1] = H_[1,5]
                H_[2,3] = EV_[0,0]*EV_[1,0]*EV_[4,0]*EV_[5,0]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[5,0]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]
                H_[5,2] = H_[2,5]
                H_[3,4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[5,0]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]
                H_[5,3] = H_[3,5]
                H_[4,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]
                H_[5,4] = H_[4,5]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePR8(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_    = (
              EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0])
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
            g_[1] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
            g_[2] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
            g_[3] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
            g_[4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
            g_[5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]
            g_[6] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]
            g_[7] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]
            if nargout>2:
                H_ = np.zeros((8,8))
                H_[0,1] = EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[4,0] = H_[0,4]
                H_[0,5] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]
                H_[5,0] = H_[0,5]
                H_[0,6] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]
                H_[6,0] = H_[0,6]
                H_[0,7] = EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]
                H_[7,0] = H_[0,7]
                H_[1,2] = EV_[0,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[4,1] = H_[1,4]
                H_[1,5] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]
                H_[5,1] = H_[1,5]
                H_[1,6] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]
                H_[6,1] = H_[1,6]
                H_[1,7] = EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]
                H_[7,1] = H_[1,7]
                H_[2,3] = EV_[0,0]*EV_[1,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]
                H_[5,2] = H_[2,5]
                H_[2,6] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]
                H_[6,2] = H_[2,6]
                H_[2,7] = EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]
                H_[7,2] = H_[2,7]
                H_[3,4] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]
                H_[5,3] = H_[3,5]
                H_[3,6] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]
                H_[6,3] = H_[3,6]
                H_[3,7] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]
                H_[7,3] = H_[3,7]
                H_[4,5] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[6,0]*EV_[7,0]
                H_[5,4] = H_[4,5]
                H_[4,6] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[7,0]
                H_[6,4] = H_[4,6]
                H_[4,7] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]
                H_[7,4] = H_[4,7]
                H_[5,6] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[7,0]
                H_[6,5] = H_[5,6]
                H_[5,7] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]
                H_[7,5] = H_[5,7]
                H_[6,7] = EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]
                H_[7,6] = H_[6,7]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePR10(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_    = (
              EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0]  = (
                  EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[1]  = (
                  EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[2]  = (
                  EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[3]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[4]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[5]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[6]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
            g_[7]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
            g_[8]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
            g_[9]  = (
                  EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,1]  = (
                      EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[2,0] = H_[0,2]
                H_[0,3]  = (
                      EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[3,0] = H_[0,3]
                H_[0,4]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[4,0] = H_[0,4]
                H_[0,5]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[5,0] = H_[0,5]
                H_[0,6]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,0] = H_[0,6]
                H_[0,7]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,0] = H_[0,7]
                H_[0,8]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,0] = H_[0,8]
                H_[0,9]  = (
                      EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,0] = H_[0,9]
                H_[1,2]  = (
                      EV_[0,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[2,1] = H_[1,2]
                H_[1,3]  = (
                      EV_[0,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[3,1] = H_[1,3]
                H_[1,4]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[4,1] = H_[1,4]
                H_[1,5]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[5,1] = H_[1,5]
                H_[1,6]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,1] = H_[1,6]
                H_[1,7]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,1] = H_[1,7]
                H_[1,8]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,1] = H_[1,8]
                H_[1,9]  = (
                      EV_[0,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,1] = H_[1,9]
                H_[2,3]  = (
                      EV_[0,0]*EV_[1,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[3,2] = H_[2,3]
                H_[2,4]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[4,2] = H_[2,4]
                H_[2,5]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[5,2] = H_[2,5]
                H_[2,6]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,2] = H_[2,6]
                H_[2,7]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,2] = H_[2,7]
                H_[2,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,2] = H_[2,8]
                H_[2,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,2] = H_[2,9]
                H_[3,4]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[4,3] = H_[3,4]
                H_[3,5]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[5,3] = H_[3,5]
                H_[3,6]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,3] = H_[3,6]
                H_[3,7]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,3] = H_[3,7]
                H_[3,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,3] = H_[3,8]
                H_[3,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,3] = H_[3,9]
                H_[4,5]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[6,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[5,4] = H_[4,5]
                H_[4,6]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,4] = H_[4,6]
                H_[4,7]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,4] = H_[4,7]
                H_[4,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,4] = H_[4,8]
                H_[4,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[5,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,4] = H_[4,9]
                H_[5,6]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[7,0]*EV_[8,0]*EV_[9,0])
                H_[6,5] = H_[5,6]
                H_[5,7]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[8,0]*EV_[9,0])
                H_[7,5] = H_[5,7]
                H_[5,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[9,0])
                H_[8,5] = H_[5,8]
                H_[5,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[6,0]*EV_[7,0]*EV_[8,0])
                H_[9,5] = H_[5,9]
                H_[6,7]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[8,0]*EV_[9,0])
                H_[7,6] = H_[6,7]
                H_[6,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[9,0])
                H_[8,6] = H_[6,8]
                H_[6,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[7,0]*EV_[8,0])
                H_[9,6] = H_[6,9]
                H_[7,8]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[9,0])
                H_[8,7] = H_[7,8]
                H_[7,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[8,0])
                H_[9,7] = H_[7,9]
                H_[8,9]  = (
                      EV_[0,0]*EV_[1,0]*EV_[2,0]*EV_[3,0]*EV_[4,0]*EV_[5,0]*EV_[6,0]*EV_[7,0])
                H_[9,8] = H_[8,9]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

