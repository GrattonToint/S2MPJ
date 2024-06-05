from s2mpjlib import *
class  DITTERT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DITTERT
#    *********
# 
#    Minimize the Dittert function.
#    Source: See Minc, Linear and Multilinear Algebra 21, 1987
# 
#    SIF input: N. Gould, March 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "OQR2-AN-V-V"
# 
#    Size of matrix
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   3              $-PARAMETER
# IE N                   4              $-PARAMETER
# IE N                   5              $-PARAMETER
# IE N                   6              $-PARAMETER
# IE N                   7              $-PARAMETER
# IE N                   8              $-PARAMETER     original value
# IE N                   9              $-PARAMETER
# IE N                   10             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DITTERT'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DITTERT'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['RN'] = float(v_['N'])
        v_['N+1'] = 1+v_['N']
        v_['2**N'] = 1
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['N-I+1'] = v_['N+1']-I
            v_['I-1'] = - 1+I
            v_['R2**N'] = float(v_['2**N'])
            v_['S'+str(int(v_['N-I+1']))] = 0.1+v_['R2**N']
            v_['T'+str(int(v_['I-1']))] = 0.1+v_['R2**N']
            v_['2**N'] = v_['2**N']*v_['2']
        v_['N-1'] = - 1+v_['N']
        v_['2**N-1'] = - 1+v_['2**N']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                [iv,ix_,_] = s2mpj_ii('P'+str(K),ix_)
                pb.xnames=arrset(pb.xnames,iv,'P'+str(K))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('A'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'A'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('R'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'R'+str(I))
            [iv,ix_,_] = s2mpj_ii('C'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'C'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('ROWPROD',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('COLPROD',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['P'+str(int(v_['2**N-1']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                [ig,ig_,_] = s2mpj_ii('PE'+str(K),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'PE'+str(K))
                iv = ix_['P'+str(K)]
                pbm.A[ig,iv] = float(- 1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('R'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'R'+str(I))
                iv = ix_['A'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('C'+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'C'+str(J))
                iv = ix_['A'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('R'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(I))
            iv = ix_['R'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['C'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('SUM',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'SUM')
            iv = ix_['R'+str(I)]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['SUM'],float(v_['RN']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['R'+str(I)]] = 1.0e-6
            pb.xlower[ix_['C'+str(I)]] = 1.0e-6
            for J in range(int(v_['1']),int(v_['N'])+1):
                pb.xupper[ix_['A'+str(I)+','+str(J)]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'P')
        [it,iet_,_] = s2mpj_ii( 'eLOG', iet_)
        elftv = loaset(elftv,it,0,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                v_['ID'] = 0
                v_['PT'] = 1
                v_['KK'] = K
                for I in range(int(v_['1']),int(v_['N'])+1):
                    v_['SI'] = v_['S'+str(I)]
                    v_['ISI'] = int(np.fix(v_['SI']))
                    v_['BI'] = int(np.fix(v_['KK']/v_['ISI']))
                    v_['ID'] = v_['ID']+v_['BI']
                    v_['BISI'] = v_['BI']*v_['ISI']
                    v_['KK'] = v_['KK']-v_['BISI']
                    v_['RI'] = float(I)
                    v_['RNZ'+str(int(v_['PT']))] = 0.1+v_['RI']
                    v_['PT'] = v_['PT']+v_['BI']
                v_['I1'] = v_['0']
                v_['I2'] = v_['1']
                v_['ID-2'] = - 2+v_['ID']
                for I in range(int(v_['1']),int(v_['ID-2'])+1):
                    v_['I1'] = v_['ID']
                    v_['I2'] = v_['0']
                for I in range(int(v_['1']),int(v_['I1'])+1):
                    v_['RJ'] = v_['RNZ'+str(I)]
                    v_['J'] = int(np.fix(v_['RJ']))
                    v_['SI'] = v_['S'+str(int(v_['J']))]
                    v_['ISI'] = int(np.fix(v_['SI']))
                    v_['IPP'] = K-v_['ISI']
                    ename = 'E'+str(K)+','+str(I)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    vname = 'A'+str(int(v_['ID']))+','+str(int(v_['J']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='A')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'P'+str(int(v_['IPP']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='P')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                for I in range(int(v_['1']),int(v_['I2'])+1):
                    v_['RJ'] = v_['RNZ'+str(int(v_['1']))]
                    v_['J'] = int(np.fix(v_['RJ']))
                    v_['RJJ'] = v_['RNZ'+str(int(v_['2']))]
                    v_['JJ'] = int(np.fix(v_['RJJ']))
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['2']))+','+str(int(v_['J']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='A')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['1']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['1']))+','+str(int(v_['JJ']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='P')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['2']))+','+str(int(v_['JJ']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='A')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    ename = 'E'+str(K)+','+str(int(v_['2']))
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    vname = 'A'+str(int(v_['1']))+','+str(int(v_['J']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='P')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                v_['RD'] = float(v_['ID'])
                v_['D'+str(K)] = 0.1+v_['RD']
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'LOGC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eLOG')
            ielftype = arrset(ielftype, ie, iet_["eLOG"])
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'LOGR'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eLOG')
            ielftype = arrset(ielftype, ie, iet_["eLOG"])
            vname = 'R'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gEXP',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for M in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RK1'] = v_['T'+str(M)]
            v_['K1'] = int(np.fix(v_['RK1']))
            v_['K2'] = 2*v_['K1']
            v_['K1'] = 1+v_['K1']
            v_['K2'] = - 1+v_['K2']
            for K in range(int(v_['K1']),int(v_['K2'])+1):
                v_['RD'] = v_['D'+str(K)]
                v_['ID'] = int(np.fix(v_['RD']))
                for I in range(int(v_['1']),int(v_['ID'])+1):
                    ig = ig_['PE'+str(K)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(K)+','+str(I)])
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['ROWPROD']
            pbm.grftype = arrset(pbm.grftype,ig,'gEXP')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['LOGR'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['COLPROD']
            pbm.grftype = arrset(pbm.grftype,ig,'gEXP')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['LOGC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(2)           5.0D-1
# LO SOLTN(3)           2.22222222D-1
# LO SOLTN(4)           9.375-2
# LO SOLTN(5)           3.84D-2
# LO SOLTN(6)           1.54321098D-2
# LO SOLTN(7)           6.11989902D-3
# LO SOLTN(8)           2.40325927D-3
# LO SOLTN(9)           9.36656708D-4
# LO SOLTN(10)          3.6288D-4
# LO SOLTN(11)          1.39905948D-4
# LO SOLTN(12)          5.37232170D-5
# LO SOLTN(13)          2.05596982D-5
# LO SOLTN(14)          7.84541375D-6
# LO SOLTN(15)          2.98628137D-6
# LO SOLTN(16)          1.13422671D-6
# LO SOLTN(17)          4.29968709D-7
# LO SOLTN(18)          1.62718123D-7
# LO SOLTN(19)          6.14859946D-8
# LO SOLTN(20)          2.32019615D-8
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
        pb.pbclass = "OQR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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
    def eLOG(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.log(EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -1.0/(EV_[0]*EV_[0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gEXP(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        EXZ = -np.exp(GVAR_)
        f_= EXZ
        if nargout>1:
            g_ = EXZ
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = EXZ
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

