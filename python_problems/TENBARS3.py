from s2mpjlib import *
class  TENBARS3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TENBARS3
#    *********
#    The "ten bar truss" structural optimization problem,
#    version P3.
# 
#    The problem is to minimize the cross section areas of the bars
#    in the structure
# 
#      /|
#      /|>o------o------o
#      /|  \    /|\    /|
#           \  / | \  / |
#            \/  |  \/  |
#            /\  |  /\  |
#           /  \ | /  \ |
#      /|  /    \|/    \|
#      /|>o------o------o
#      /|
# 
#    submitted to vertical forces of equal magnitude (P0) applied at
#    the two free lower nodes, subject to limits of nodal displacements.
# 
#    Source:
#    K. Svanberg,
#    private communication,  August 1990.
#    See also
#    K. Svanberg,
#    "On local and global minima in structural optimization",
#    in "New directions in optimum structural design" (Atrek, Ragsdell
#    and Zienkiwewicz, eds.), Wiley, 1984.
# 
#    SIF input: Ph. Toint, August 1990.
# 
#    classification = "C-CLOR2-MY-18-8"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TENBARS3'

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
        v_['8'] = 8
        v_['10'] = 10
        v_['SQ2'] = np.sqrt(2.0)
        v_['SQ8'] = np.sqrt(8.0)
        v_['1/SQ8'] = 1.0/v_['SQ8']
        v_['-1/SQ8'] = -1.0*v_['1/SQ8']
        v_['C0'] = 2.53106
        v_['-P0'] = -589.884
        v_['C0SQ2'] = v_['C0']*v_['SQ2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['8'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        for I in range(int(v_['1']),int(v_['10'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X1']])
        valA = np.append(valA,float(v_['C0']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X2']])
        valA = np.append(valA,float(v_['C0SQ2']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X3']])
        valA = np.append(valA,float(v_['C0SQ2']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X4']])
        valA = np.append(valA,float(v_['C0']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X5']])
        valA = np.append(valA,float(v_['C0']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X6']])
        valA = np.append(valA,float(v_['C0']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X7']])
        valA = np.append(valA,float(v_['C0SQ2']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X8']])
        valA = np.append(valA,float(v_['C0SQ2']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X9']])
        valA = np.append(valA,float(v_['C0']))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X10']])
        valA = np.append(valA,float(v_['C0']))
        for I in range(int(v_['1']),int(v_['8'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        self.gconst = arrset(self.gconst,ig_['C4'],float(v_['-P0']))
        self.gconst = arrset(self.gconst,ig_['C8'],float(v_['-P0']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['10'])+1):
            self.xlower[ix_['X'+str(I)]] = 0.645
        self.xlower[ix_['U4']] = -50.8
        self.xlower[ix_['U8']] = -50.8
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXU', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'U')
        [it,iet_,_] = s2mpj_ii( 'eXUPV', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'U')
        elftv = loaset(elftv,it,2,'V')
        [it,iet_,_] = s2mpj_ii( 'eXUMV', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'U')
        elftv = loaset(elftv,it,2,'V')
        [it,iet_,_] = s2mpj_ii( 'eXBIG', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'U')
        elftv = loaset(elftv,it,2,'V')
        elftv = loaset(elftv,it,3,'W')
        elftv = loaset(elftv,it,4,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'EA'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXU')
        ielftype = arrset(ielftype,ie,iet_["eXU"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EB'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXUPV')
        ielftype = arrset(ielftype,ie,iet_["eXUPV"])
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EC'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eXUMV')
            ielftype = arrset(ielftype,ie,iet_['eXUMV'])
        vname = 'X6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'ED'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXBIG')
        ielftype = arrset(ielftype,ie,iet_["eXBIG"])
        vname = 'X7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EE'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eXUMV')
            ielftype = arrset(ielftype,ie,iet_['eXUMV'])
        vname = 'X5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EF'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eXUMV')
            ielftype = arrset(ielftype,ie,iet_['eXUMV'])
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EG'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXU')
        ielftype = arrset(ielftype,ie,iet_["eXU"])
        vname = 'X4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EH'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eXBIG')
        ielftype = arrset(ielftype,ie,iet_["eXBIG"])
        vname = 'X8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='Z')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EI'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eXUMV')
            ielftype = arrset(ielftype,ie,iet_['eXUMV'])
        vname = 'X9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'EJ'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eXUMV')
            ielftype = arrset(ielftype,ie,iet_['eXUMV'])
        vname = 'X10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='U')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'U8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(0.0))
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['C1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EA'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EC'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EB'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['ED'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        ig = ig_['C2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EB'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EE'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['ED'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/SQ8']))
        ig = ig_['C3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EF'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EI'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['EG'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['C4']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EF'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EE'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C5']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EC'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C6']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EJ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['C7']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['ED'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['-1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EI'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        ig = ig_['C8']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['ED'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['1/SQ8']))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['EJ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               2247.1290
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
        self.pbclass   = "C-CLOR2-MY-18-8"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXU(self, nargout,*args):

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
    def eXUPV(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        U_[1,0] = U_[1,0]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        f_   = IV_[1]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[1] = IV_[0]
            g_[0] = IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = 1.0
                H_[0,1] = H_[1,0]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXUMV(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]-1
        U_[1,0] = U_[1,0]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        f_   = IV_[1]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[1] = IV_[0]
            g_[0] = IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = 1.0
                H_[0,1] = H_[1,0]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXBIG(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,5))
        IV_ = np.zeros(2)
        U_[0,1] = U_[0,1]+1
        U_[0,2] = U_[0,2]+1
        U_[0,3] = U_[0,3]-1
        U_[0,4] = U_[0,4]-1
        U_[1,0] = U_[1,0]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        f_   = IV_[1]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[1] = IV_[0]
            g_[0] = IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[1,0] = 1.0
                H_[0,1] = H_[1,0]
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

