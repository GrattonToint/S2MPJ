from s2mpjlib import *
class  WATSONNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : WATSONNE
#    *********
# 
#    Watson problem in 12 variables. This is a nonlinear equation version
#    of problem WATSON.
# 
#    This function  is a nonlinear least squares with 31 groups.  Each
#    group has 1 nonlinear and 1 linear elements.
# 
#    Source:  problem 20 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#128 (p. 100).
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "NOR2-AN-V-31"
# 
#    The number of variables can be varied, but should be smaller than
#    31
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   12             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'WATSONNE'

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
        if nargin<1:
            v_['N'] = int(12);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   31             $-PARAMETER
        v_['M'] = 31
        v_['1'] = 1
        v_['2'] = 2
        v_['29'] = 29
        v_['30'] = 30
        v_['1/29'] = 0.024482759
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['29'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['1/29']
            v_['LNTI'] = np.log(v_['TI'])
            for J in range(int(v_['2']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['RJ-1'] = -1.0+v_['RJ']
                v_['RJ-2'] = -2.0+v_['RJ']
                v_['AE'] = v_['RJ-2']*v_['LNTI']
                v_['C0'] = np.exp(v_['AE'])
                v_['C'] = v_['C0']*v_['RJ-1']
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I))
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['30'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G'+str(int(v_['30'])))
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G'+str(int(v_['M'])))
        iv = ix_['X'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['G'+str(int(v_['30']))],float(0.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eMWSQ', iet_)
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
        elftv = loaset(elftv,it,10,'V11')
        elftv = loaset(elftv,it,11,'V12')
        elftp = []
        elftp = loaset(elftp,it,0,'T1')
        elftp = loaset(elftp,it,1,'T2')
        elftp = loaset(elftp,it,2,'T3')
        elftp = loaset(elftp,it,3,'T4')
        elftp = loaset(elftp,it,4,'T5')
        elftp = loaset(elftp,it,5,'T6')
        elftp = loaset(elftp,it,6,'T7')
        elftp = loaset(elftp,it,7,'T8')
        elftp = loaset(elftp,it,8,'T9')
        elftp = loaset(elftp,it,9,'T10')
        elftp = loaset(elftp,it,10,'T11')
        elftp = loaset(elftp,it,11,'T12')
        [it,iet_,_] = s2mpj_ii( 'eMSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['29'])+1):
            v_['RI'] = float(I)
            v_['TI'] = v_['RI']*v_['1/29']
            v_['LNTI'] = np.log(v_['TI'])
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eMWSQ')
            ielftype = arrset(ielftype, ie, iet_["eMWSQ"])
            pb.x0 = np.zeros((pb.n,1))
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
            vname = 'X11'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V11')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X12'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V12')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['1']),int(v_['N'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                v_['CE0'] = v_['RJ-1']*v_['LNTI']
                v_['CE'+str(J)] = np.exp(v_['CE0'])
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T1')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE1']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T2')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE2']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T3')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE3']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T4')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE4']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T5')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE5']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T6')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE6']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T7')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE7']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T8')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE8']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T9')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE9']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T10')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE10']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T11')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE11']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='T12')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['CE12']))
        ename = 'E'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eMSQ')
        ielftype = arrset(ielftype, ie, iet_["eMSQ"])
        ename = 'E'+str(int(v_['M']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['29'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G'+str(int(v_['M']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN(12)           2.27559922D-9
# LO SOLTN(31)           1.53795068D-9
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
        pb.pbclass = "NOR2-AN-V-31"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eMSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = -EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EV_[0]-EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eMWSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U  = (
              pbm.elpar[iel_][0]*EV_[0]+pbm.elpar[iel_][1]*EV_[1]+pbm.elpar[iel_][2]*EV_[2]+pbm.elpar[iel_][3]*EV_[3]+pbm.elpar[iel_][4]*EV_[4]+pbm.elpar[iel_][5]*EV_[5]+pbm.elpar[iel_][6]*EV_[6]+pbm.elpar[iel_][7]*EV_[7]+pbm.elpar[iel_][8]*EV_[8]+pbm.elpar[iel_][9]*EV_[9]+pbm.elpar[iel_][10]*EV_[10]+pbm.elpar[iel_][11]*EV_[11])
        TWOT1 = pbm.elpar[iel_][0]+pbm.elpar[iel_][0]
        TWOT2 = pbm.elpar[iel_][1]+pbm.elpar[iel_][1]
        TWOT3 = pbm.elpar[iel_][2]+pbm.elpar[iel_][2]
        TWOT4 = pbm.elpar[iel_][3]+pbm.elpar[iel_][3]
        TWOT5 = pbm.elpar[iel_][4]+pbm.elpar[iel_][4]
        TWOT6 = pbm.elpar[iel_][5]+pbm.elpar[iel_][5]
        TWOT7 = pbm.elpar[iel_][6]+pbm.elpar[iel_][6]
        TWOT8 = pbm.elpar[iel_][7]+pbm.elpar[iel_][7]
        TWOT9 = pbm.elpar[iel_][8]+pbm.elpar[iel_][8]
        TWOT10 = pbm.elpar[iel_][9]+pbm.elpar[iel_][9]
        TWOT11 = pbm.elpar[iel_][10]+pbm.elpar[iel_][10]
        TWOT12 = pbm.elpar[iel_][11]+pbm.elpar[iel_][11]
        f_   = -U*U
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -TWOT1*U
            g_[1] = -TWOT2*U
            g_[2] = -TWOT3*U
            g_[3] = -TWOT4*U
            g_[4] = -TWOT5*U
            g_[5] = -TWOT6*U
            g_[6] = -TWOT7*U
            g_[7] = -TWOT8*U
            g_[8] = -TWOT9*U
            g_[9] = -TWOT10*U
            g_[10] = -TWOT11*U
            g_[11] = -TWOT12*U
            if nargout>2:
                H_ = np.zeros((12,12))
                H_[0,0] = -TWOT1*pbm.elpar[iel_][0]
                H_[0,1] = -TWOT1*pbm.elpar[iel_][1]
                H_[1,0] = H_[0,1]
                H_[0,2] = -TWOT1*pbm.elpar[iel_][2]
                H_[2,0] = H_[0,2]
                H_[0,3] = -TWOT1*pbm.elpar[iel_][3]
                H_[3,0] = H_[0,3]
                H_[0,4] = -TWOT1*pbm.elpar[iel_][4]
                H_[4,0] = H_[0,4]
                H_[0,5] = -TWOT1*pbm.elpar[iel_][5]
                H_[5,0] = H_[0,5]
                H_[0,6] = -TWOT1*pbm.elpar[iel_][6]
                H_[6,0] = H_[0,6]
                H_[0,7] = -TWOT1*pbm.elpar[iel_][7]
                H_[7,0] = H_[0,7]
                H_[0,8] = -TWOT1*pbm.elpar[iel_][8]
                H_[8,0] = H_[0,8]
                H_[0,9] = -TWOT1*pbm.elpar[iel_][9]
                H_[9,0] = H_[0,9]
                H_[0,10] = -TWOT1*pbm.elpar[iel_][10]
                H_[10,0] = H_[0,10]
                H_[0,11] = -TWOT1*pbm.elpar[iel_][11]
                H_[11,0] = H_[0,11]
                H_[1,1] = -TWOT2*pbm.elpar[iel_][1]
                H_[1,2] = -TWOT2*pbm.elpar[iel_][2]
                H_[2,1] = H_[1,2]
                H_[1,3] = -TWOT2*pbm.elpar[iel_][3]
                H_[3,1] = H_[1,3]
                H_[1,4] = -TWOT2*pbm.elpar[iel_][4]
                H_[4,1] = H_[1,4]
                H_[1,5] = -TWOT2*pbm.elpar[iel_][5]
                H_[5,1] = H_[1,5]
                H_[1,6] = -TWOT2*pbm.elpar[iel_][6]
                H_[6,1] = H_[1,6]
                H_[1,7] = -TWOT2*pbm.elpar[iel_][7]
                H_[7,1] = H_[1,7]
                H_[1,8] = -TWOT2*pbm.elpar[iel_][7]
                H_[8,1] = H_[1,8]
                H_[1,9] = -TWOT2*pbm.elpar[iel_][9]
                H_[9,1] = H_[1,9]
                H_[1,10] = -TWOT2*pbm.elpar[iel_][10]
                H_[10,1] = H_[1,10]
                H_[1,11] = -TWOT2*pbm.elpar[iel_][11]
                H_[11,1] = H_[1,11]
                H_[2,2] = -TWOT3*pbm.elpar[iel_][2]
                H_[2,3] = -TWOT3*pbm.elpar[iel_][3]
                H_[3,2] = H_[2,3]
                H_[2,4] = -TWOT3*pbm.elpar[iel_][4]
                H_[4,2] = H_[2,4]
                H_[2,5] = -TWOT3*pbm.elpar[iel_][5]
                H_[5,2] = H_[2,5]
                H_[2,6] = -TWOT3*pbm.elpar[iel_][6]
                H_[6,2] = H_[2,6]
                H_[2,7] = -TWOT3*pbm.elpar[iel_][7]
                H_[7,2] = H_[2,7]
                H_[2,8] = -TWOT3*pbm.elpar[iel_][7]
                H_[8,2] = H_[2,8]
                H_[2,9] = -TWOT3*pbm.elpar[iel_][9]
                H_[9,2] = H_[2,9]
                H_[2,10] = -TWOT3*pbm.elpar[iel_][10]
                H_[10,2] = H_[2,10]
                H_[2,11] = -TWOT3*pbm.elpar[iel_][11]
                H_[11,2] = H_[2,11]
                H_[3,3] = -TWOT4*pbm.elpar[iel_][3]
                H_[3,4] = -TWOT4*pbm.elpar[iel_][4]
                H_[4,3] = H_[3,4]
                H_[3,5] = -TWOT4*pbm.elpar[iel_][5]
                H_[5,3] = H_[3,5]
                H_[3,6] = -TWOT4*pbm.elpar[iel_][6]
                H_[6,3] = H_[3,6]
                H_[3,7] = -TWOT4*pbm.elpar[iel_][7]
                H_[7,3] = H_[3,7]
                H_[3,8] = -TWOT4*pbm.elpar[iel_][7]
                H_[8,3] = H_[3,8]
                H_[3,9] = -TWOT4*pbm.elpar[iel_][9]
                H_[9,3] = H_[3,9]
                H_[3,10] = -TWOT4*pbm.elpar[iel_][10]
                H_[10,3] = H_[3,10]
                H_[3,11] = -TWOT4*pbm.elpar[iel_][11]
                H_[11,3] = H_[3,11]
                H_[4,4] = -TWOT5*pbm.elpar[iel_][4]
                H_[4,5] = -TWOT5*pbm.elpar[iel_][5]
                H_[5,4] = H_[4,5]
                H_[4,6] = -TWOT5*pbm.elpar[iel_][6]
                H_[6,4] = H_[4,6]
                H_[4,7] = -TWOT5*pbm.elpar[iel_][7]
                H_[7,4] = H_[4,7]
                H_[4,8] = -TWOT5*pbm.elpar[iel_][7]
                H_[8,4] = H_[4,8]
                H_[4,9] = -TWOT5*pbm.elpar[iel_][9]
                H_[9,4] = H_[4,9]
                H_[4,10] = -TWOT5*pbm.elpar[iel_][10]
                H_[10,4] = H_[4,10]
                H_[4,11] = -TWOT5*pbm.elpar[iel_][11]
                H_[11,4] = H_[4,11]
                H_[5,5] = -TWOT6*pbm.elpar[iel_][5]
                H_[5,6] = -TWOT6*pbm.elpar[iel_][6]
                H_[6,5] = H_[5,6]
                H_[5,7] = -TWOT6*pbm.elpar[iel_][7]
                H_[7,5] = H_[5,7]
                H_[5,8] = -TWOT6*pbm.elpar[iel_][7]
                H_[8,5] = H_[5,8]
                H_[5,9] = -TWOT6*pbm.elpar[iel_][9]
                H_[9,5] = H_[5,9]
                H_[5,10] = -TWOT6*pbm.elpar[iel_][10]
                H_[10,5] = H_[5,10]
                H_[5,11] = -TWOT6*pbm.elpar[iel_][11]
                H_[11,5] = H_[5,11]
                H_[6,6] = -TWOT7*pbm.elpar[iel_][6]
                H_[6,7] = -TWOT7*pbm.elpar[iel_][7]
                H_[7,6] = H_[6,7]
                H_[6,8] = -TWOT7*pbm.elpar[iel_][7]
                H_[8,6] = H_[6,8]
                H_[6,9] = -TWOT7*pbm.elpar[iel_][9]
                H_[9,6] = H_[6,9]
                H_[6,10] = -TWOT7*pbm.elpar[iel_][10]
                H_[10,6] = H_[6,10]
                H_[6,11] = -TWOT7*pbm.elpar[iel_][11]
                H_[11,6] = H_[6,11]
                H_[7,7] = -TWOT8*pbm.elpar[iel_][7]
                H_[7,8] = -TWOT8*pbm.elpar[iel_][7]
                H_[8,7] = H_[7,8]
                H_[7,9] = -TWOT8*pbm.elpar[iel_][9]
                H_[9,7] = H_[7,9]
                H_[7,10] = -TWOT8*pbm.elpar[iel_][10]
                H_[10,7] = H_[7,10]
                H_[7,11] = -TWOT8*pbm.elpar[iel_][11]
                H_[11,7] = H_[7,11]
                H_[8,8] = -TWOT9*pbm.elpar[iel_][8]
                H_[8,9] = -TWOT9*pbm.elpar[iel_][9]
                H_[9,8] = H_[8,9]
                H_[8,10] = -TWOT9*pbm.elpar[iel_][10]
                H_[10,8] = H_[8,10]
                H_[8,11] = -TWOT9*pbm.elpar[iel_][11]
                H_[11,8] = H_[8,11]
                H_[9,9] = -TWOT10*pbm.elpar[iel_][9]
                H_[9,10] = -TWOT10*pbm.elpar[iel_][10]
                H_[10,9] = H_[9,10]
                H_[9,11] = -TWOT10*pbm.elpar[iel_][11]
                H_[11,9] = H_[9,11]
                H_[10,10] = -TWOT11*pbm.elpar[iel_][10]
                H_[10,11] = -TWOT11*pbm.elpar[iel_][11]
                H_[11,10] = H_[10,11]
                H_[11,11] = -TWOT12*pbm.elpar[iel_][11]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

