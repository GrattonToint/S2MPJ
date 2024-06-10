from s2mpjlib import *
class  ROSZMAN1LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ROSZMAN1LS
#    *********
# 
#    NIST Data fitting problem ROSZMAN1.
# 
#    Fit: y =  b1 - b2*x - arctan[b3/(x-b4)]/pi + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#   Reference: Roszman, L., NIST (1979).  
#     Quantum Defects for Sulfur I Atom.
# 
#    classification = "SUR2-MN-4-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ROSZMAN1LS'

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
        v_['M'] = 25
        v_['N'] = 4
        v_['1'] = 1
        v_['X1'] = -4868.68
        v_['X2'] = -4868.09
        v_['X3'] = -4867.41
        v_['X4'] = -3375.19
        v_['X5'] = -3373.14
        v_['X6'] = -3372.03
        v_['X7'] = -2473.74
        v_['X8'] = -2472.35
        v_['X9'] = -2469.45
        v_['X10'] = -1894.65
        v_['X11'] = -1893.40
        v_['X12'] = -1497.24
        v_['X13'] = -1495.85
        v_['X14'] = -1493.41
        v_['X15'] = -1208.68
        v_['X16'] = -1206.18
        v_['X17'] = -1206.04
        v_['X18'] = -997.92
        v_['X19'] = -996.61
        v_['X20'] = -996.31
        v_['X21'] = -834.94
        v_['X22'] = -834.66
        v_['X23'] = -710.03
        v_['X24'] = -530.16
        v_['X25'] = -464.17
        v_['Y1'] = 0.252429
        v_['Y2'] = 0.252141
        v_['Y3'] = 0.251809
        v_['Y4'] = 0.297989
        v_['Y5'] = 0.296257
        v_['Y6'] = 0.295319
        v_['Y7'] = 0.339603
        v_['Y8'] = 0.337731
        v_['Y9'] = 0.333820
        v_['Y10'] = 0.389510
        v_['Y11'] = 0.386998
        v_['Y12'] = 0.438864
        v_['Y13'] = 0.434887
        v_['Y14'] = 0.427893
        v_['Y15'] = 0.471568
        v_['Y16'] = 0.461699
        v_['Y17'] = 0.461144
        v_['Y18'] = 0.513532
        v_['Y19'] = 0.506641
        v_['Y20'] = 0.505062
        v_['Y21'] = 0.535648
        v_['Y22'] = 0.533726
        v_['Y23'] = 0.568064
        v_['Y24'] = 0.612886
        v_['Y25'] = 0.624169
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['B1']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            v_['-X'] = -1.0*v_['X'+str(I)]
            iv = ix_['B2']
            pbm.A[ig,iv] = float(v_['-X'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['B1']] = float(0.1)
        pb.x0[ix_['B2']] = float(-0.00001)
        pb.x0[ix_['B3']] = float(1000.0)
        pb.x0[ix_['B4']] = float(-100.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE7', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE7')
            ielftype = arrset(ielftype, ie, iet_["eE7"])
            vname = 'B3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['X'+str(I)]))
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
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-MN-4-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,4.0*np.arctan(1.0e0))
        return pbm

    @staticmethod
    def eE7(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V12 = EV_[0]*EV_[0]
        V13 = EV_[0]*V12
        V2MX = EV_[1]-pbm.elpar[iel_][0]
        V2MX2 = V2MX*V2MX
        V2MX3 = V2MX*V2MX2
        R = V12/V2MX2+1.0
        PIR = pbm.efpar[0]*R
        PIR2 = PIR*R
        f_   = -np.arctan(EV_[0]/V2MX)/pbm.efpar[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/(PIR*V2MX)
            g_[1] = EV_[0]/(PIR*V2MX2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*EV_[0]/(PIR2*V2MX3)
                H_[0,1] = 1.0/(PIR*V2MX2)-2.0*V12/(PIR2*V2MX**4)
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*V13/(PIR2*V2MX**5)-2.0*EV_[0]/(PIR*V2MX3)
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

