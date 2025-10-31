from s2mpjlib import *
class  SPMSRTLS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPMSRTLS
#    *********
# 
#    Liu and Nocedal tridiagonal matrix square root problem.
# 
#    Source:  problem 151 (p. 93) in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    This is a least-squares variant of problem SPMSQRT.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CSUR2-AN-V-V"
# 
#    M is the dimension of the matrix
#    The number of variables is 3*M-2
# 
#           Alternative values for the SIF file parameters:
# IE M                   10             $-PARAMETER n = 28     original value
# IE M                   34             $-PARAMETER n = 100
# IE M                   167            $-PARAMETER n = 499
# IE M                   334            $-PARAMETER n = 1000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPMSRTLS'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(1667);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
# IE M                   3334           $-PARAMETER n = 10000
        v_['M-1'] = -1+v_['M']
        v_['M-2'] = -2+v_['M']
        v_['M-3'] = -3+v_['M']
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['B1,1'] = np.sin(1.0)
        v_['B1,2'] = np.sin(4.0)
        v_['K'] = 2
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['K'] = 1+v_['K']
            v_['RK'] = float(v_['K'])
            v_['RKSQR'] = v_['RK']*v_['RK']
            v_['B'+str(I)+','+str(int(v_['I-1']))] = np.sin(v_['RKSQR'])
            v_['K'] = 1+v_['K']
            v_['RK'] = float(v_['K'])
            v_['RKSQR'] = v_['RK']*v_['RK']
            v_['B'+str(I)+','+str(I)] = np.sin(v_['RKSQR'])
            v_['K'] = 1+v_['K']
            v_['RK'] = float(v_['K'])
            v_['RKSQR'] = v_['RK']*v_['RK']
            v_['B'+str(I)+','+str(int(v_['I+1']))] = np.sin(v_['RKSQR'])
        v_['K'] = 1+v_['K']
        v_['RK'] = float(v_['K'])
        v_['RKSQR'] = v_['RK']*v_['RK']
        v_['B'+str(int(v_['M']))+','+str(int(v_['M-1']))] = np.sin(v_['RKSQR'])
        v_['K'] = 1+v_['K']
        v_['RK'] = float(v_['K'])
        v_['RKSQR'] = v_['RK']*v_['RK']
        v_['B'+str(int(v_['M']))+','+str(int(v_['M']))] = np.sin(v_['RKSQR'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('X'+str(int(v_['1']))+','+str(int(v_['1'])),ix_)
        self.xnames = (
             arrset(self.xnames,iv,'X'+str(int(v_['1']))+','+str(int(v_['1']))))
        [iv,ix_,_] = s2mpj_ii('X'+str(int(v_['1']))+','+str(int(v_['2'])),ix_)
        self.xnames = (
             arrset(self.xnames,iv,'X'+str(int(v_['1']))+','+str(int(v_['2']))))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(int(v_['I-1'])),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(int(v_['I-1'])))
            [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(I))
            [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(int(v_['I+1'])),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(int(v_['I+1'])))
        [iv,ix_,_] = s2mpj_ii('X'+str(int(v_['M']))+','+str(int(v_['M-1'])),ix_)
        self.xnames = (
             arrset(self.xnames,iv,'X'+str(int(v_['M']))+','+str(int(v_['M-1']))))
        [iv,ix_,_] = s2mpj_ii('X'+str(int(v_['M']))+','+str(int(v_['M'])),ix_)
        self.xnames = (
             arrset(self.xnames,iv,'X'+str(int(v_['M']))+','+str(int(v_['M']))))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for J in range(int(v_['1']),int(v_['3'])+1):
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['1']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
        for J in range(int(v_['1']),int(v_['4'])+1):
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['2']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['3']),int(v_['M-2'])+1):
            v_['I-2'] = -2+I
            v_['I+2'] = 2+I
            for J in range(int(v_['I-2']),int(v_['I+2'])+1):
                [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        for J in range(int(v_['M-3']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['M-1']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
        for J in range(int(v_['M-2']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['M']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        v_['ENTRY'] = v_['B1,1']*v_['B1,1']
        v_['PROD'] = v_['B1,2']*v_['B2,1']
        v_['ENTRY'] = v_['ENTRY']+v_['PROD']
        self.gconst = arrset(self.gconst,ig_['E1,1'],float(v_['ENTRY']))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['ENTRY'] = v_['B'+str(I)+','+str(I)]*v_['B'+str(I)+','+str(I)]
            v_['PROD'] = (v_['B'+str(int(v_['I-1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(int(v_['I-1']))])
            v_['ENTRY'] = v_['ENTRY']+v_['PROD']
            v_['PROD'] = (v_['B'+str(int(v_['I+1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(int(v_['I+1']))])
            v_['ENTRY'] = v_['ENTRY']+v_['PROD']
            self.gconst  = (
                  arrset(self.gconst,ig_['E'+str(I)+','+str(I)],float(v_['ENTRY'])))
        v_['ENTRY'] = (v_['B'+str(int(v_['M']))+','+str(int(v_['M']))]*v_['B'+
             str(int(v_['M']))+','+str(int(v_['M']))])
        v_['PROD'] = (v_['B'+str(int(v_['M-1']))+','+str(int(v_['M']))]*v_['B'+
             str(int(v_['M']))+','+str(int(v_['M-1']))])
        v_['ENTRY'] = v_['ENTRY']+v_['PROD']
        self.gconst  = (
              arrset(self.gconst,ig_['E'+str(int(v_['M']))+','+str(int(v_['M']))],float(v_['ENTRY'])))
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            v_['I+1'] = 1+I
            v_['ENTRY'] = (v_['B'+str(int(v_['I+1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(I)])
            v_['PROD'] = (v_['B'+str(int(v_['I+1']))+','+str(int(v_['I+1']))]*v_['B'+
                 str(int(v_['I+1']))+','+str(I)])
            v_['ENTRY'] = v_['ENTRY']+v_['PROD']
            self.gconst  = (
                  arrset(self.gconst,ig_['E'+str(int(v_['I+1']))+','+str(I)],float(v_['ENTRY'])))
        for I in range(int(v_['2']),int(v_['M'])+1):
            v_['I-1'] = -1+I
            v_['ENTRY'] = (v_['B'+str(int(v_['I-1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(I)])
            v_['PROD'] = (v_['B'+str(int(v_['I-1']))+','+str(int(v_['I-1']))]*v_['B'+
                 str(int(v_['I-1']))+','+str(I)])
            v_['ENTRY'] = v_['ENTRY']+v_['PROD']
            self.gconst  = (
                  arrset(self.gconst,ig_['E'+str(int(v_['I-1']))+','+str(I)],float(v_['ENTRY'])))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['ENTRY'] = (v_['B'+str(int(v_['I+1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(int(v_['I-1']))])
            self.gconst  = (
                  arrset(self.gconst,ig_['E'+str(int(v_['I+1']))+','+str(int(v_['I-1']))],float(v_['ENTRY'])))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['ENTRY'] = (v_['B'+str(int(v_['I-1']))+','+str(I)]*v_['B'+str(I)+
                 ','+str(int(v_['I+1']))])
            self.gconst  = (
                  arrset(self.gconst,ig_['E'+str(int(v_['I-1']))+','+str(int(v_['I+1']))],float(v_['ENTRY'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        v_['PROD'] = 0.2*v_['B1,1']
        self.x0[ix_['X1,1']] = float(v_['PROD'])
        v_['PROD'] = 0.2*v_['B1,2']
        self.x0[ix_['X1,2']] = float(v_['PROD'])
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['PROD'] = 0.2*v_['B'+str(I)+','+str(int(v_['I-1']))]
            self.x0[ix_['X'+str(I)+','+str(int(v_['I-1']))]] = float(v_['PROD'])
            v_['PROD'] = 0.2*v_['B'+str(I)+','+str(I)]
            self.x0[ix_['X'+str(I)+','+str(I)]] = float(v_['PROD'])
            v_['PROD'] = 0.2*v_['B'+str(I)+','+str(int(v_['I+1']))]
            self.x0[ix_['X'+str(I)+','+str(int(v_['I+1']))]] = float(v_['PROD'])
        v_['PROD'] = 0.2*v_['B'+str(int(v_['M']))+','+str(int(v_['M-1']))]
        self.x0[ix_['X'+str(int(v_['M']))+','+str(int(v_['M-1']))]]  = (
              float(v_['PROD']))
        v_['PROD'] = 0.2*v_['B'+str(int(v_['M']))+','+str(int(v_['M']))]
        self.x0[ix_['X'+str(int(v_['M']))+','+str(int(v_['M']))]]  = (
              float(v_['PROD']))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD2', iet_)
        elftv = loaset(elftv,it,0,'VI')
        elftv = loaset(elftv,it,1,'VJ')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'D'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        ename = 'D'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'G'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'H'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'R'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'S'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'B'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'C'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'F'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        ename = 'D'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'G'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['3']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'H'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'R'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['3']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'S'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S'+str(int(v_['2']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['3']))+','+str(int(v_['4']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['3']),int(v_['M-2'])+1):
            v_['I-2'] = -2+I
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))+','+str(int(v_['I-2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))+','+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'F'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I-1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'D'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'X'+str(I)+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'G'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'H'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)+','+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'R'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD2')
            ielftype = arrset(ielftype,ie,iet_["ePROD2"])
            vname = 'X'+str(I)+','+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VI')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))+','+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'A'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-2']))+','+str(int(v_['M-3']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'B'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-2']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'C'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'F'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-2']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        ename = 'D'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'G'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'G'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'H'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'H'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'R'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'R'+str(int(v_['M-1']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'A'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'A'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-2']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'B'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'B'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'C'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'ePROD2')
        ielftype = arrset(ielftype,ie,iet_["ePROD2"])
        ename = 'F'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M-1']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VI')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'F'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M-1']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='VJ')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSQ')
        ielftype = arrset(ielftype,ie,iet_["eSQ"])
        ename = 'D'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X'+str(int(v_['M']))+','+str(int(v_['M']))
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
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
        ig = ig_['E'+str(int(v_['1']))+','+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['G'+str(int(v_['1']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['1']))+','+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['H'+str(int(v_['1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['R'+str(int(v_['1']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['1']))+','+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['E'+str(int(v_['2']))+','+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['2']))+','+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['G'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['2']))+','+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['H'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['R'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['2']))+','+str(int(v_['4']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(int(v_['2']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['3']),int(v_['M-2'])+1):
            v_['I-2'] = -2+I
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            ig = ig_['E'+str(I)+','+str(int(v_['I-2']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['E'+str(I)+','+str(int(v_['I-1']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['E'+str(I)+','+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['F'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['G'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['E'+str(I)+','+str(int(v_['I+1']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['H'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['R'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel, 1.)
            ig = ig_['E'+str(I)+','+str(int(v_['I+2']))]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['E'+str(int(v_['M-1']))+','+str(int(v_['M-3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['E'+str(int(v_['M-1']))+','+str(int(v_['M-2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['M-1']))+','+str(int(v_['M-1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['G'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['M-1']))+','+str(int(v_['M']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['H'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['R'+str(int(v_['M-1']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['M']))+','+str(int(v_['M-2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(int(v_['M']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['E'+str(int(v_['M']))+','+str(int(v_['M-1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['B'+str(int(v_['M']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C'+str(int(v_['M']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        ig = ig_['E'+str(int(v_['M']))+','+str(int(v_['M']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['F'+str(int(v_['M']))])
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['D'+str(int(v_['M']))])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CSUR2-AN-V-V"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD2(self, nargout,*args):

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
    def eSQ(self, nargout,*args):

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

