from s2mpjlib import *
class  LUBRIF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUBRIF
#    *********
# 
#    The elastodynamic lubrification problem by Kostreva.
#    NB. This version has an error. 
#    See LUBRIFC for the correct formulation.
# 
#    Source:
#    M.M. Kostreva,
#    "Elasto-hydrodynamic lubrification: a non-linear
#    complementarity problem",
#    International Journal for Numerical Methods in Fluids,
#    4: 377-397, 1984.
# 
#    This problem is problem #5 in More's test set.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "QOR2-MN-V-V"
# 
#    Number of discretized points per unit length
# 
#           Alternative values for the SIF file parameters:
# IE NN                  10             $-PARAMETER n = 151    original value
# IE NN                  50             $-PARAMETER n = 751
# IE NN                  250            $-PARAMETER n = 3751
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LUBRIF'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LUBRIF'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NN'] = int(5);  #  SIF file default value
        else:
            v_['NN'] = int(args[0])
        v_['ALPHA'] = 1.838
        v_['LAMBDA'] = 1.642
        v_['XA'] = -3.0
        v_['XF'] = 2.0
        v_['N'] = 5*v_['NN']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['PI'] = 3.1415926535
        v_['2N'] = 2*v_['N']
        v_['2N-2'] = -2+v_['2N']
        v_['2N-1'] = -1+v_['2N']
        v_['2N+2'] = 2+v_['2N']
        v_['-XA'] = -1.0*v_['XA']
        v_['LEN'] = v_['XF']+v_['-XA']
        v_['1/PI'] = 1.0/v_['PI']
        v_['1/2PI'] = 0.5*v_['1/PI']
        v_['RN'] = float(v_['N'])
        v_['1/N'] = 1.0/v_['RN']
        v_['DX'] = v_['LEN']*v_['1/N']
        v_['1/DX'] = 1.0/v_['DX']
        v_['L/DX'] = v_['LAMBDA']*v_['1/DX']
        v_['-L/DX'] = -1.0*v_['L/DX']
        v_['1/DX2'] = v_['1/DX']*v_['1/DX']
        v_['-1/DX2'] = -1.0*v_['1/DX2']
        v_['DX/PI'] = v_['DX']*v_['1/PI']
        v_['2DX/PI'] = 2.0*v_['DX/PI']
        v_['DX/2'] = 0.5*v_['DX']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('K',ix_)
        pb.xnames=arrset(pb.xnames,iv,'K')
        for I in range(int(v_['0']),int(v_['2N'])+1,int(v_['2'])):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
        for J in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            [iv,ix_,_] = s2mpj_ii('H'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'H'+str(J))
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            [iv,ix_,_] = s2mpj_ii('R'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'R'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(int(v_['0'])))
            iv = ix_['P'+str(I)]
            pbm.A[ig,iv] = float(v_['2DX/PI'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'R'+str(int(v_['0'])))
        iv = ix_['P'+str(int(v_['2N']))]
        pbm.A[ig,iv] = float(v_['DX/PI'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('COMPL',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('DR'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'DR'+str(I))
            iv = ix_['H'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['L/DX'])+pbm.A[ig,iv]
            iv = ix_['H'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-L/DX'])+pbm.A[ig,iv]
            iv = ix_['R'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            v_['-J'] = -1*J
            [ig,ig_,_] = s2mpj_ii('DH'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'DH'+str(J))
            iv = ix_['K']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['H'+str(J)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for I in range(int(v_['2']),int(v_['2N'])+1):
                v_['C'+str(I)] = 0.0
            v_['RI-J'] = float(v_['-J'])
            v_['I-JDX'] = v_['RI-J']*v_['DX/2']
            v_['ALN'] = np.absolute(v_['I-JDX'])
            v_['LN'] = np.log(v_['ALN'])
            v_['T1'] = v_['I-JDX']*v_['LN']
            v_['COEFF'] = v_['T1']*v_['1/2PI']
            v_['C'+str(int(v_['2']))] = v_['C'+str(int(v_['2']))]+v_['COEFF']
            v_['I-J'] = 2+v_['-J']
            v_['RI-J'] = float(v_['I-J'])
            v_['I-JDX'] = v_['RI-J']*v_['DX/2']
            v_['ALN'] = np.absolute(v_['I-JDX'])
            v_['LN'] = np.log(v_['ALN'])
            v_['T1'] = v_['I-JDX']*v_['LN']
            v_['COEFF'] = v_['T1']*v_['1/PI']
            v_['C'+str(int(v_['4']))] = v_['C'+str(int(v_['4']))]+v_['COEFF']
            for I in range(int(v_['4']),int(v_['2N-2'])+1,int(v_['2'])):
                v_['I-2'] = -2+I
                v_['I+2'] = 2+I
                v_['I-J'] = I+v_['-J']
                v_['RI-J'] = float(v_['I-J'])
                v_['I-JDX'] = v_['RI-J']*v_['DX/2']
                v_['ALN'] = np.absolute(v_['I-JDX'])
                v_['LN'] = np.log(v_['ALN'])
                v_['T1'] = v_['I-JDX']*v_['LN']
                v_['COEFF'] = v_['T1']*v_['1/PI']
                v_['C'+str(int(v_['I+2']))] = v_['C'+str(int(v_['I+2']))]+v_['COEFF']
                v_['-COEFF'] = -1.0*v_['COEFF']
                v_['C'+str(int(v_['I-2']))] = v_['C'+str(int(v_['I-2']))]+v_['-COEFF']
            v_['I-J'] = v_['2N']+v_['-J']
            v_['RI-J'] = float(v_['I-J'])
            v_['I-JDX'] = v_['RI-J']*v_['DX/2']
            v_['ALN'] = np.absolute(v_['I-JDX'])
            v_['LN'] = np.log(v_['ALN'])
            v_['T1'] = v_['I-JDX']*v_['LN']
            v_['COEFF'] = v_['T1']*v_['1/2PI']
            v_['-COEFF'] = -1.0*v_['COEFF']
            v_['C'+str(int(v_['2N-2']))] = v_['C'+str(int(v_['2N-2']))]+v_['-COEFF']
            for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
                [ig,ig_,_] = s2mpj_ii('DH'+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'DH'+str(J))
                iv = ix_['P'+str(I)]
                pbm.A[ig,iv] = float(v_['C'+str(I)])+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['R'+str(int(v_['0']))],float(1.0))
        for J in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            v_['RJ'] = float(J)
            v_['JDX'] = v_['RJ']*v_['DX/2']
            v_['XJ'] = v_['XA']+v_['JDX']
            v_['XJSQ'] = v_['XJ']*v_['XJ']
            v_['XJSQ+1'] = 1.0+v_['XJSQ']
            v_['RHS'] = -1.0*v_['XJSQ+1']
            pbm.gconst = arrset(pbm.gconst,ig_['DH'+str(J)],float(v_['RHS']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['K']] = -float('Inf')
        pb.xupper[ix_['K']] = +float('Inf')
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            pb.xupper[ix_['P'+str(I)]] = 3.0
            pb.xlower[ix_['P'+str(I)]] = 0.0
        for I in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            pb.xlower[ix_['H'+str(I)]] = -float('Inf')
            pb.xupper[ix_['H'+str(I)]] = +float('Inf')
        pb.xlower[ix_['P'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['P'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['P'+str(int(v_['2N']))]] = 0.0
        pb.xupper[ix_['P'+str(int(v_['2N']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        v_['2NN'] = v_['NN']+v_['NN']
        v_['4NN'] = 4*v_['NN']
        for I in range(int(v_['2']),int(v_['4NN'])+1,int(v_['2'])):
            v_['RI'] = float(I)
            v_['IDX'] = v_['RI']*v_['DX/2']
            v_['XI'] = v_['XA']+v_['IDX']
            v_['LIN'] = 0.02*v_['XI']
            v_['PI0'] = 0.06+v_['LIN']
            if('P'+str(I) in ix_):
                pb.x0[ix_['P'+str(I)]] = float(v_['PI0'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(I)]),float(v_['PI0'])))
        v_['4NN+2'] = 2+v_['4NN']
        v_['8NN'] = 8*v_['NN']
        for I in range(int(v_['4NN+2']),int(v_['8NN'])+1,int(v_['2'])):
            v_['RI'] = float(I)
            v_['IDX'] = v_['RI']*v_['DX/2']
            v_['XI'] = v_['XA']+v_['IDX']
            v_['XISQ'] = v_['XI']*v_['XI']
            v_['-XISQ'] = -1.0*v_['XISQ']
            v_['1-XISQ'] = 1.0+v_['-XISQ']
            v_['PI0'] = np.sqrt(v_['1-XISQ'])
            if('P'+str(I) in ix_):
                pb.x0[ix_['P'+str(I)]] = float(v_['PI0'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(I)]),float(v_['PI0'])))
        v_['8NN+2'] = 2+v_['8NN']
        for I in range(int(v_['8NN+2']),int(v_['2N'])+1,int(v_['2'])):
            if('P'+str(I) in ix_):
                pb.x0[ix_['P'+str(I)]] = float(0.0)
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P'+str(I)]),float(0.0)))
        for J in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            v_['RJ'] = float(J)
            v_['JDX'] = v_['RJ']*v_['DX/2']
            v_['XJ'] = v_['XA']+v_['JDX']
            v_['XJSQ'] = v_['XJ']*v_['XJ']
            if('H'+str(J) in ix_):
                pb.x0[ix_['H'+str(J)]] = float(v_['XJSQ'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['H'+str(J)]),float(v_['XJSQ'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eREY', iet_)
        elftv = loaset(elftv,it,0,'PA')
        elftv = loaset(elftv,it,1,'PB')
        elftv = loaset(elftv,it,2,'H')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'P')
        elftv = loaset(elftv,it,1,'R')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for J in range(int(v_['1']),int(v_['2N-1'])+1,int(v_['2'])):
            v_['I+'] = 1+J
            v_['I-'] = -1+J
            ename = 'ER'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eREY')
            ielftype = arrset(ielftype, ie, iet_["eREY"])
            vname = 'P'+str(int(v_['I-']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='PA')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(int(v_['I+']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='H')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'H'+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='PB')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['ALPHA']))
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'R'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='R')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            ig = ig_['COMPL']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['2N-2'])+1,int(v_['2'])):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ig = ig_['DR'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ER'+str(int(v_['I-1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/DX2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ER'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/DX2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
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
        pb.pbclass = "QOR2-MN-V-V"
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
    def eREY(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        HA = -0.5*pbm.elpar[iel_][0]
        EARG = HA*(EV_[0]+EV_[1])
        E = np.exp(EARG)
        PAMPB = EV_[0]-EV_[1]
        T1 = PAMPB*HA+1.0
        T2 = PAMPB*HA-1.0
        HSQ = EV_[2]*EV_[2]
        HCB = HSQ*EV_[2]
        f_   = PAMPB*HCB*E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = T1*HCB*E
            g_[1] = T2*HCB*E
            g_[2] = 3.0*PAMPB*HSQ*E
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = HCB*E*HA*(T1+1.0)
                H_[0,1] = HCB*E*HA*(T1-1.0)
                H_[1,0] = H_[0,1]
                H_[0,2] = 3.0*T1*HSQ*E
                H_[2,0] = H_[0,2]
                H_[1,1] = HCB*E*HA*(T2-1.0)
                H_[1,2] = 3.0*T2*HSQ*E
                H_[2,1] = H_[1,2]
                H_[2,2] = 6.0*EV_[2]*PAMPB*E
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

