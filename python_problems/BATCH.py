from s2mpjlib import *
class  BATCH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BATCH
#    *********
# 
#    Source: Optimal Design of Multiproduct Batch Plant
#    G.R. Kocis & I.E. Grossmann,
#    "Global OPtimization of Nonconvex Mixed Integer Nonlinear Programmming
#     (MINLP) problems in Process Synthesis", Indust. Engng. Chem. Res.,
#    No. 27, pp 1407--1421, 1988.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "OOR2-AN-46-73"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BATCH'

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
        v_['M'] = 6
        v_['N'] = 5
        v_['NU'] = 4
        v_['LOGNU'] = np.log(4.0)
        v_['VL'] = np.log(300.0)
        v_['VU'] = np.log(3000.0)
        v_['H'] = 6000.0
        v_['TLO1'] = 0.729961
        v_['TLO2'] = 0.530628
        v_['TLO3'] = 1.09024
        v_['TLO4'] = -0.133531
        v_['TLO5'] = 0.0487901
        v_['TUP1'] = 2.11626
        v_['TUP2'] = 1.91626
        v_['TUP3'] = 2.47654
        v_['TUP4'] = 1.25276
        v_['TUP5'] = 1.43508
        v_['BLO1'] = 4.45966
        v_['BLO2'] = 3.74950
        v_['BLO3'] = 4.49144
        v_['BLO4'] = 3.14988
        v_['BLO5'] = 3.04452
        v_['BUP1'] = 397.747
        v_['BUP2'] = 882.353
        v_['BUP3'] = 833.333
        v_['BUP4'] = 638.298
        v_['BUP5'] = 666.667
        v_['Q1'] = 250000.0
        v_['Q2'] = 150000.0
        v_['Q3'] = 180000.0
        v_['Q4'] = 160000.0
        v_['Q5'] = 120000.0
        v_['LOGI1'] = np.log(1.0)
        v_['LOGI2'] = np.log(2.0)
        v_['LOGI3'] = np.log(3.0)
        v_['LOGI4'] = np.log(4.0)
        v_['S1,1'] = np.log(7.9)
        v_['S2,1'] = np.log(0.7)
        v_['S3,1'] = np.log(0.7)
        v_['S4,1'] = np.log(4.7)
        v_['S5,1'] = np.log(1.2)
        v_['S1,2'] = np.log(2.0)
        v_['S2,2'] = np.log(0.8)
        v_['S3,2'] = np.log(2.6)
        v_['S4,2'] = np.log(2.3)
        v_['S5,2'] = np.log(3.6)
        v_['S1,3'] = np.log(5.2)
        v_['S2,3'] = np.log(0.9)
        v_['S3,3'] = np.log(1.6)
        v_['S4,3'] = np.log(1.6)
        v_['S5,3'] = np.log(2.4)
        v_['S1,4'] = np.log(4.9)
        v_['S2,4'] = np.log(3.4)
        v_['S3,4'] = np.log(3.6)
        v_['S4,4'] = np.log(2.7)
        v_['S5,4'] = np.log(4.5)
        v_['S1,5'] = np.log(6.1)
        v_['S2,5'] = np.log(2.1)
        v_['S3,5'] = np.log(3.2)
        v_['S4,5'] = np.log(1.2)
        v_['S5,5'] = np.log(1.6)
        v_['S1,6'] = np.log(4.2)
        v_['S2,6'] = np.log(2.5)
        v_['S3,6'] = np.log(2.9)
        v_['S4,6'] = np.log(2.5)
        v_['S5,6'] = np.log(2.1)
        v_['T1,1'] = np.log(6.4)
        v_['T2,1'] = np.log(6.8)
        v_['T3,1'] = np.log(1.0)
        v_['T4,1'] = np.log(3.2)
        v_['T5,1'] = np.log(2.1)
        v_['T1,2'] = np.log(4.7)
        v_['T2,2'] = np.log(6.4)
        v_['T3,2'] = np.log(6.3)
        v_['T4,2'] = np.log(3.0)
        v_['T5,2'] = np.log(2.5)
        v_['T1,3'] = np.log(8.3)
        v_['T2,3'] = np.log(6.5)
        v_['T3,3'] = np.log(5.4)
        v_['T4,3'] = np.log(3.5)
        v_['T5,3'] = np.log(4.2)
        v_['T1,4'] = np.log(3.9)
        v_['T2,4'] = np.log(4.4)
        v_['T3,4'] = np.log(11.9)
        v_['T4,4'] = np.log(3.3)
        v_['T5,4'] = np.log(3.6)
        v_['T1,5'] = np.log(2.1)
        v_['T2,5'] = np.log(2.3)
        v_['T3,5'] = np.log(5.7)
        v_['T4,5'] = np.log(2.8)
        v_['T5,5'] = np.log(3.7)
        v_['T1,6'] = np.log(1.2)
        v_['T2,6'] = np.log(3.2)
        v_['T3,6'] = np.log(6.2)
        v_['T4,6'] = np.log(3.4)
        v_['T5,6'] = np.log(2.2)
        for J in range(int(v_['1']),int(v_['M'])+1):
            v_['ALPHA'+str(J)] = 250.0
            v_['BETA'+str(J)] = 0.6
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            [iv,ix_,_] = s2mpj_ii('N'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'N'+str(J))
        for J in range(int(v_['1']),int(v_['M'])+1):
            [iv,ix_,_] = s2mpj_ii('V'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(J))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('TL'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'TL'+str(I))
        for J in range(int(v_['1']),int(v_['M'])+1):
            for K in range(int(v_['1']),int(v_['NU'])+1):
                [iv,ix_,_] = s2mpj_ii('Y'+str(K)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(K)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('COST',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('VOL'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'VOL'+str(I)+','+str(J))
                iv = ix_['V'+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['B'+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('CYCL'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CYCL'+str(I)+','+str(J))
                iv = ix_['N'+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['TL'+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('HORIZON',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'HORIZON')
        for J in range(int(v_['1']),int(v_['M'])+1):
            for K in range(int(v_['1']),int(v_['NU'])+1):
                [ig,ig_,_] = s2mpj_ii('NPAR'+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'NPAR'+str(J))
                iv = ix_['Y'+str(K)+','+str(J)]
                pbm.A[ig,iv] = float(v_['LOGI'+str(K)])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('NPAR'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'NPAR'+str(J))
            iv = ix_['N'+str(J)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['M'])+1):
            for K in range(int(v_['1']),int(v_['NU'])+1):
                [ig,ig_,_] = s2mpj_ii('SOS1'+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'SOS1'+str(J))
                iv = ix_['Y'+str(K)+','+str(J)]
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['VOL'+str(I)+','+str(J)],float(v_['S'+str(I)+','+str(J)])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['CYCL'+str(I)+','+str(J)],float(v_['T'+str(I)+','+str(J)])))
        pbm.gconst = arrset(pbm.gconst,ig_['HORIZON'],float(v_['H']))
        for J in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['SOS1'+str(J)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for J in range(int(v_['1']),int(v_['M'])+1):
            pb.xupper[ix_['N'+str(J)]] = v_['LOGNU']
            pb.xlower[ix_['V'+str(J)]] = v_['VL']
            pb.xupper[ix_['V'+str(J)]] = v_['VU']
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['B'+str(I)]] = v_['BLO'+str(I)]
            pb.xupper[ix_['B'+str(I)]] = v_['BUP'+str(I)]
            pb.xlower[ix_['TL'+str(I)]] = v_['TLO'+str(I)]
            pb.xupper[ix_['TL'+str(I)]] = v_['TUP'+str(I)]
        for J in range(int(v_['1']),int(v_['M'])+1):
            for K in range(int(v_['1']),int(v_['NU'])+1):
                pb.xupper[ix_['Y'+str(K)+','+str(J)]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXPXAY', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for J in range(int(v_['1']),int(v_['M'])+1):
            ename = 'EXPO'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXPXAY')
            ielftype = arrset(ielftype, ie, iet_["eEXPXAY"])
            pb.x0 = np.zeros((pb.n,1))
            vname = 'N'+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'V'+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['BETA'+str(J)]))
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'EXPC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXPXAY')
            ielftype = arrset(ielftype, ie, iet_["eEXPXAY"])
            vname = 'TL'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(-1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['COST']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPO'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['ALPHA'+str(J)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['HORIZON']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EXPC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['Q'+str(I)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-46-73"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXPXAY(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        FVALUE = np.exp(EV_[0]+pbm.elpar[iel_][0]*EV_[1])
        GYVALU = pbm.elpar[iel_][0]*FVALUE
        f_   = FVALUE
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = FVALUE
            g_[1] = GYVALU
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = FVALUE
                H_[0,1] = GYVALU
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.elpar[iel_][0]*GYVALU
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

