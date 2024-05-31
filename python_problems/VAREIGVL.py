from s2xlib import *
class  VAREIGVL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : VAREIGVL
#    *********
# 
#    The variational eigenvalue by Auchmuty.
#    This problems features a banded matrix of bandwidth 2M+1 = 9.
# 
#    This problem has N least-squares groups, each having a linear part
#    only and N nonlinear elements,
#    plus a least q-th power group having N nonlinear elements.
# 
#    Source: problem 1 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995
#               and Nick Gould, December, 2019, May 2024
# 
#    classification = "OUR2-AN-V-0"
# 
#    Number of variables -1 (variable)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'VAREIGVL'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'VAREIGVL'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(19);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   49             $-PARAMETER     original value
# IE N                   99             $-PARAMETER
# IE N                   499            $-PARAMETER
# IE N                   999            $-PARAMETER
# IE N                   4999           $-PARAMETER
# IE M                   4              $-PARAMETER  .le. N
# IE M                   5              $-PARAMETER  .le. N
        if nargin<2:
            v_['M'] = int(6);  #  SIF file default value
        else:
            v_['M'] = int(args[1])
        if nargin<3:
            v_['Q'] = float(1.5);  #  SIF file default value
        else:
            v_['Q'] = float(args[2])
        v_['1'] = 1
        v_['-1.0'] = -1.0
        v_['N+1'] = 1+v_['N']
        v_['-M'] = -1*v_['M']
        v_['M+1'] = 1+v_['M']
        v_['N-M'] = v_['N']+v_['-M']
        v_['N-M+1'] = 1+v_['N-M']
        v_['N2'] = v_['N']*v_['N']
        v_['RN2'] = float(v_['N2'])
        v_['-1/N2'] = v_['-1.0']/v_['RN2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        [iv,ix_,_] = s2x_ii('MU',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MU')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I+M'] = I+v_['M']
            for J in range(int(v_['1']),int(v_['I+M'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['AIJ'])+pbm.A[ig,iv]
        for I in range(int(v_['M+1']),int(v_['N-M'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I-M'] = I+v_['-M']
            v_['I+M'] = I+v_['M']
            for J in range(int(v_['I-M']),int(v_['I+M'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['AIJ'])+pbm.A[ig,iv]
        for I in range(int(v_['N-M+1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['-I'] = -1.0*v_['RI']
            v_['I-M'] = I+v_['-M']
            for J in range(int(v_['I-M']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['IJ'] = v_['RI']*v_['RJ']
                v_['SIJ'] = np.sin(v_['IJ'])
                v_['J-I'] = v_['RJ']+v_['-I']
                v_['J-ISQ'] = v_['J-I']*v_['J-I']
                v_['ARG'] = v_['J-ISQ']*v_['-1/N2']
                v_['EX'] = np.exp(v_['ARG'])
                v_['AIJ'] = v_['SIJ']*v_['EX']
                [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['AIJ'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('G'+str(int(v_['N+1'])),ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        pb.x0[ix_['MU']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'M')
        elftv = loaset(elftv,it,1,'X')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'P'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'MU'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='M')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'S'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gLQ',igt_)
        [it,igt_,_] = s2x_ii('gLQ',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'POWER')
        [it,igt_,_] = s2x_ii('gLQ2',igt_)
        [it,igt_,_] = s2x_ii('gLQ2',igt_)
        grftp = loaset(grftp,it,0,'POWER')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gLQ')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='POWER')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(2.0))
        ig = ig_['G'+str(int(v_['N+1']))]
        pbm.grftype = arrset(pbm.grftype,ig,'gLQ2')
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['G'+str(int(v_['N+1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['G'+str(int(v_['N+1']))]
        posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='POWER')
        pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['Q']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(pbm,nargout,*args):

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gLQ(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        IPOWER = pbm.grpar[igr_][0]
        PM1 = IPOWER-1
        f_= GVAR_**IPOWER/pbm.grpar[igr_][0]
        if nargout>1:
            g_ = GVAR_**PM1
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = PM1*GVAR_**(IPOWER-2)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gLQ2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_**pbm.grpar[igr_][0]/pbm.grpar[igr_][0]
        if nargout>1:
            g_ = GVAR_**(pbm.grpar[igr_][0]-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = (pbm.grpar[igr_][0]-1.0e0)*GVAR_**(pbm.grpar[igr_][0]-2.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

