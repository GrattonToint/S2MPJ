from s2mpjlib import *
class  AIRPORT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This problem is concerned with the localisation of airports in Brazil.
#    We consider  m  balls in the real plane, whose centers are the coordinates
#    of some Brazilian  cities and whose  radius were chosen such that the balls are
#    disjoint. The problem is to find one point  (xi, yi) on  each ball, i=1,..,m,
#    such that  SUM(||(xi,yi) - (xj,yj)||)  is  minimum, where the sum involves all
#    the pairs (i,j) such that 1 <= i <= m, 1 <= j <= m and i <> j.
# 
#    For this problem instance, we have m =  42 cities and n = 84 points, 
#    i.e, 42 nonlinear inequalities constraints and 84 variables.
# 
#    Source:
#    Contribution from a LANCELOT user.
# 
#    SIF input : Rodrigo de Barros Nabholz & Maria Aparecida Diniz Ehrhardt
#                November 1994, DMA - IMECC- UNICAMP
#    Adaptation for CUTE: Ph. Toint, November 1994.
# 
#    classification = "SQR2-MN-84-42"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'AIRPORT'

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
        v_['N'] = 42
        v_['N-1'] = 41
        v_['1'] = 1
        v_['R1'] = 0.09
        v_['R2'] = 0.3
        v_['R3'] = 0.09
        v_['R4'] = 0.45
        v_['R5'] = 0.5
        v_['R6'] = 0.04
        v_['R7'] = 0.1
        v_['R8'] = 0.02
        v_['R9'] = 0.02
        v_['R10'] = 0.07
        v_['R11'] = 0.4
        v_['R12'] = 0.045
        v_['R13'] = 0.05
        v_['R14'] = 0.056
        v_['R15'] = 0.36
        v_['R16'] = 0.08
        v_['R17'] = 0.07
        v_['R18'] = 0.36
        v_['R19'] = 0.67
        v_['R20'] = 0.38
        v_['R21'] = 0.37
        v_['R22'] = 0.05
        v_['R23'] = 0.4
        v_['R24'] = 0.66
        v_['R25'] = 0.05
        v_['R26'] = 0.07
        v_['R27'] = 0.08
        v_['R28'] = 0.3
        v_['R29'] = 0.31
        v_['R30'] = 0.49
        v_['R31'] = 0.09
        v_['R32'] = 0.46
        v_['R33'] = 0.12
        v_['R34'] = 0.07
        v_['R35'] = 0.07
        v_['R36'] = 0.09
        v_['R37'] = 0.05
        v_['R38'] = 0.13
        v_['R39'] = 0.16
        v_['R40'] = 0.46
        v_['R41'] = 0.25
        v_['R42'] = 0.1
        v_['CX1'] = -6.3
        v_['CX2'] = -7.8
        v_['CX3'] = -9.0
        v_['CX4'] = -7.2
        v_['CX5'] = -5.7
        v_['CX6'] = -1.9
        v_['CX7'] = -3.5
        v_['CX8'] = -0.5
        v_['CX9'] = 1.4
        v_['CX10'] = 4.0
        v_['CX11'] = 2.1
        v_['CX12'] = 5.5
        v_['CX13'] = 5.7
        v_['CX14'] = 5.7
        v_['CX15'] = 3.8
        v_['CX16'] = 5.3
        v_['CX17'] = 4.7
        v_['CX18'] = 3.3
        v_['CX19'] = 0.0
        v_['CX20'] = -1.0
        v_['CX21'] = -0.4
        v_['CX22'] = 4.2
        v_['CX23'] = 3.2
        v_['CX24'] = 1.7
        v_['CX25'] = 3.3
        v_['CX26'] = 2.0
        v_['CX27'] = 0.7
        v_['CX28'] = 0.1
        v_['CX29'] = -0.1
        v_['CX30'] = -3.5
        v_['CX31'] = -4.0
        v_['CX32'] = -2.7
        v_['CX33'] = -0.5
        v_['CX34'] = -2.9
        v_['CX35'] = -1.2
        v_['CX36'] = -0.4
        v_['CX37'] = -0.1
        v_['CX38'] = -1.0
        v_['CX39'] = -1.7
        v_['CX40'] = -2.1
        v_['CX41'] = -1.8
        v_['CX42'] = 0.0
        v_['CY1'] = 8.0
        v_['CY2'] = 5.1
        v_['CY3'] = 2.0
        v_['CY4'] = 2.6
        v_['CY5'] = 5.5
        v_['CY6'] = 7.1
        v_['CY7'] = 5.9
        v_['CY8'] = 6.6
        v_['CY9'] = 6.1
        v_['CY10'] = 5.6
        v_['CY11'] = 4.9
        v_['CY12'] = 4.7
        v_['CY13'] = 4.3
        v_['CY14'] = 3.6
        v_['CY15'] = 4.1
        v_['CY16'] = 3.0
        v_['CY17'] = 2.4
        v_['CY18'] = 3.0
        v_['CY19'] = 4.7
        v_['CY20'] = 3.4
        v_['CY21'] = 2.3
        v_['CY22'] = 1.5
        v_['CY23'] = 0.5
        v_['CY24'] = -1.7
        v_['CY25'] = -2.0
        v_['CY26'] = -3.1
        v_['CY27'] = -3.5
        v_['CY28'] = -2.4
        v_['CY29'] = -1.3
        v_['CY30'] = 0.0
        v_['CY31'] = -1.7
        v_['CY32'] = -2.1
        v_['CY33'] = -0.4
        v_['CY34'] = -2.9
        v_['CY35'] = -3.4
        v_['CY36'] = -4.3
        v_['CY37'] = -5.2
        v_['CY38'] = -6.5
        v_['CY39'] = -7.5
        v_['CY40'] = -6.4
        v_['CY41'] = -5.1
        v_['CY42'] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = I+v_['1']
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ1'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('OBJ2'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Y'+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['Y'+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('CONS'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CONS'+str(I))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CONS'+str(I)],float(v_['R'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['X'+str(I)]] = -10
            pb.xupper[ix_['X'+str(I)]] = 10
            pb.xlower[ix_['Y'+str(I)]] = -10
            pb.xupper[ix_['Y'+str(I)]] = 10
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eDIFSQR', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = []
        elftp = loaset(elftp,it,0,'W')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eDIFSQR')
            ielftype = arrset(ielftype, ie, iet_["eDIFSQR"])
            pb.x0 = np.zeros((pb.n,1))
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='W')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CX'+str(I)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eDIFSQR')
            ielftype = arrset(ielftype, ie, iet_["eDIFSQR"])
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='W')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CY'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = I+v_['1']
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['OBJ1'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
                ig = ig_['OBJ2'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['CONS'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = .0
#    Solution
# LO SOLTN              47952.695811
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "SQR2-MN-84-42"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eDIFSQR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DIF = EV_[0]-pbm.elpar[iel_][0]
        f_   = DIF*DIF
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*DIF
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
    def gSQUARE(pbm,nargout,*args):

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

