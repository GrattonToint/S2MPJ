from s2mpjlib import *
class  CYCLOOCTLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The cyclooctane molecule is comprised of eight carbon atoms aligned
#    in an equally spaced ring. When they take a position of minimum
#    potential energy so that next-neighbours are equally spaced.
# 
#    Given positions v_1, ..., v_p in R^3 (with p = 8 for cyclooctane),
#    and given a spacing c^2 we have that
# 
#       ||v_i - v_i+1,mod p||^2 = c^2 for i = 1,..,p, and
#       ||v_i - v_i+2,mod p||^2 = 2p/(p-2) c^3
# 
#    where (arbitrarily) we have v_1 = 0 and component 1 of v_2 = 0
# 
#    Source:
#    an extension of the cyclooctane molecule configuration space as
#    described in (for example)
# 
#     E. Coutsias, S. Martin, A. Thompson & J. Watson
#     "Topology of cyclooctane energy landscape"
#     J. Chem. Phys. 132-234115 (2010)
# 
#    SIF input: Nick Gould, Feb 2020.
# 
#    This is a least-squares version of CYCLOOCT.SIF
# 
#    classification = "C-CSBR2-MN-V-0"
# 
#    The number of molecules
# 
#           Alternative values for the SIF file parameters:
# IE P                   8              $-PARAMETER     original value
# IE P                   100            $-PARAMETER
# IE P                   1000           $-PARAMETER
# IE P                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CYCLOOCTLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(8);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
# IE P                   100000         $-PARAMETER
        v_['C'] = 1.0
        v_['1'] = 1
        v_['2'] = 2
        v_['P-1'] = -1+v_['P']
        v_['P-2'] = -2+v_['P']
        v_['THREE'] = 3.0
        v_['RP'] = float(v_['P'])
        v_['2RP'] = 2.0*v_['RP']
        v_['RP-2'] = -2.0+v_['RP']
        v_['2RP/RP-2'] = v_['2RP']/v_['RP-2']
        v_['C2'] = v_['C']*v_['C']
        v_['SC2'] = v_['2RP/RP-2']*v_['C2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['P'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
            [iv,ix_,_] = s2mpj_ii('Z'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Z'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['P'])+1):
            [ig,ig_,_] = s2mpj_ii('A'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        selfnob      = ngrp
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['P'])+1):
            self.gconst = arrset(self.gconst,ig_['A'+str(I)],float(v_['C2']))
            self.gconst = arrset(self.gconst,ig_['B'+str(I)],float(v_['SC2']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['X1']] = 0.0
        self.xupper[ix_['X1']] = 0.0
        self.xlower[ix_['Y1']] = 0.0
        self.xupper[ix_['Y1']] = 0.0
        self.xlower[ix_['Z1']] = 0.0
        self.xupper[ix_['Z1']] = 0.0
        self.xlower[ix_['X2']] = 0.0
        self.xupper[ix_['X2']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['P'])+1):
            v_['RI'] = float(I)
            v_['START'] = v_['RI']/v_['RP']
            self.x0[ix_['X'+str(I)]] = float(v_['START'])
            self.x0[ix_['Y'+str(I)]] = float(v_['START'])
            self.x0[ix_['Z'+str(I)]] = float(v_['START'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2mpj_ii( 'eSQRDIF', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['P-1'])+1):
            v_['I+1'] = 1+I
            ename = 'AX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'AY'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'AZ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'Z'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'AX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'AY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'AZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'AZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['P-2'])+1):
            v_['I+2'] = 2+I
            ename = 'BX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'BY'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'BZ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQRDIF')
            ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
            vname = 'Z'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Z'+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BX'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BX'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['P-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BX'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BY'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BY'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['P-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BY'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BZ'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BZ'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['P-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BZ'+str(int(v_['P-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BX'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BY'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Y'+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQRDIF')
        ielftype = arrset(ielftype,ie,iet_["eSQRDIF"])
        ename = 'BZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['P']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'BZ'+str(int(v_['P']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'Z'+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V2')[0]
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
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['P'])+1):
            ig = ig_['A'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['AX'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['AY'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['AZ'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['B'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['BX'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['BY'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['BZ'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION             0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSBR2-MN-V-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(self, nargout,*args):

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
    def eSQRDIF(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        f_   = IV_[0]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
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

