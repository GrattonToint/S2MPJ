from s2xlib import *
class  CERI651C(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CERI651C
#    *********
# 
#    ISIS Data fitting problem CERI651C given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = c + l * x + I*A*B/2(A+B) *
#               [ exp( A*[A*S^2+2(x-X0)]/2) * erfc( A*S^2+(x-X0)/S*sqrt(2) ) +
#                 exp( B*[B*S^2+2(x-X0)]/2) * erfc( B*S^2+(x-X0)/S*sqrt(2) ) ]
# 
#    Source: fit to a sum of a linear background and a back-to-back exponential
#    using data enginx_ceria193749_spectrum_number_651_vana_corrected-0
#    from Mantid (http://www.mantidproject.org)
# 
#    subset X in [23919.5789114, 24189.3183142]
# 
#    SIF input: Nick Gould and Tyrone Rees, Mar 2016
# 
#    classification = "NOR2-MN-7-56"
# 
#    Potential and actual number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CERI651C'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CERI651C'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['MPOT'] = 10186
        v_['M'] = 56
        v_['MLOWER'] = 6916
        v_['MUPPER'] = 6971
        v_['N'] = 7
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['X6916'] = 23920.10938
        v_['X6917'] = 23924.89063
        v_['X6918'] = 23929.67188
        v_['X6919'] = 23934.45313
        v_['X6920'] = 23939.23438
        v_['X6921'] = 23944.01563
        v_['X6922'] = 23948.79688
        v_['X6923'] = 23953.57813
        v_['X6924'] = 23958.35938
        v_['X6925'] = 23963.14063
        v_['X6926'] = 23967.92188
        v_['X6927'] = 23972.70313
        v_['X6928'] = 23977.48438
        v_['X6929'] = 23982.26563
        v_['X6930'] = 23987.06250
        v_['X6931'] = 23991.87500
        v_['X6932'] = 23996.68750
        v_['X6933'] = 24001.50000
        v_['X6934'] = 24006.31250
        v_['X6935'] = 24011.12500
        v_['X6936'] = 24015.93750
        v_['X6937'] = 24020.75000
        v_['X6938'] = 24025.56250
        v_['X6939'] = 24030.37500
        v_['X6940'] = 24035.18750
        v_['X6941'] = 24040.00000
        v_['X6942'] = 24044.81250
        v_['X6943'] = 24049.62500
        v_['X6944'] = 24054.43750
        v_['X6945'] = 24059.25000
        v_['X6946'] = 24064.06250
        v_['X6947'] = 24068.87500
        v_['X6948'] = 24073.68750
        v_['X6949'] = 24078.50000
        v_['X6950'] = 24083.31250
        v_['X6951'] = 24088.12500
        v_['X6952'] = 24092.93750
        v_['X6953'] = 24097.75000
        v_['X6954'] = 24102.56250
        v_['X6955'] = 24107.37500
        v_['X6956'] = 24112.18750
        v_['X6957'] = 24117.00000
        v_['X6958'] = 24121.81250
        v_['X6959'] = 24126.62500
        v_['X6960'] = 24131.43750
        v_['X6961'] = 24136.25000
        v_['X6962'] = 24141.06250
        v_['X6963'] = 24145.89063
        v_['X6964'] = 24150.73438
        v_['X6965'] = 24155.57813
        v_['X6966'] = 24160.42188
        v_['X6967'] = 24165.26563
        v_['X6968'] = 24170.10938
        v_['X6969'] = 24174.95313
        v_['X6970'] = 24179.79688
        v_['X6971'] = 24184.64063
        v_['Y6916'] = 0.00000000
        v_['Y6917'] = 0.98041658
        v_['Y6918'] = 1.96083316
        v_['Y6919'] = 0.00000000
        v_['Y6920'] = 0.98041658
        v_['Y6921'] = 0.00000000
        v_['Y6922'] = 0.00000000
        v_['Y6923'] = 3.92166632
        v_['Y6924'] = 0.98041658
        v_['Y6925'] = 0.00000000
        v_['Y6926'] = 0.98041658
        v_['Y6927'] = 2.94124974
        v_['Y6928'] = 1.96083316
        v_['Y6929'] = 0.98041658
        v_['Y6930'] = 2.94124974
        v_['Y6931'] = 8.82374922
        v_['Y6932'] = 5.88249948
        v_['Y6933'] = 6.86291606
        v_['Y6934'] = 8.82374922
        v_['Y6935'] = 11.76499896
        v_['Y6936'] = 12.74541554
        v_['Y6937'] = 6.86291606
        v_['Y6938'] = 8.82374922
        v_['Y6939'] = 12.74541554
        v_['Y6940'] = 13.72583212
        v_['Y6941'] = 8.82374922
        v_['Y6942'] = 12.74541554
        v_['Y6943'] = 19.60833160
        v_['Y6944'] = 4.90208290
        v_['Y6945'] = 2.94124974
        v_['Y6946'] = 1.96083316
        v_['Y6947'] = 3.92166632
        v_['Y6948'] = 3.92166632
        v_['Y6949'] = 5.88249948
        v_['Y6950'] = 2.94124974
        v_['Y6951'] = 4.90208290
        v_['Y6952'] = 6.86291606
        v_['Y6953'] = 2.94124974
        v_['Y6954'] = 1.96083316
        v_['Y6955'] = 0.00000000
        v_['Y6956'] = 1.96083316
        v_['Y6957'] = 2.94124974
        v_['Y6958'] = 1.96083316
        v_['Y6959'] = 1.96083316
        v_['Y6960'] = 1.96083316
        v_['Y6961'] = 3.92166632
        v_['Y6962'] = 0.00000000
        v_['Y6963'] = 0.00000000
        v_['Y6964'] = 3.92166632
        v_['Y6965'] = 2.94124974
        v_['Y6966'] = 1.96083316
        v_['Y6967'] = 0.00000000
        v_['Y6968'] = 1.96083316
        v_['Y6969'] = 0.00000000
        v_['Y6970'] = 0.98041658
        v_['Y6971'] = 0.98041658
        v_['E6916'] = 1.00000000
        v_['E6917'] = 1.00000000
        v_['E6918'] = 1.41421356
        v_['E6919'] = 1.00000000
        v_['E6920'] = 1.00000000
        v_['E6921'] = 1.00000000
        v_['E6922'] = 1.00000000
        v_['E6923'] = 2.00000000
        v_['E6924'] = 1.00000000
        v_['E6925'] = 1.00000000
        v_['E6926'] = 1.00000000
        v_['E6927'] = 1.73205081
        v_['E6928'] = 1.41421356
        v_['E6929'] = 1.00000000
        v_['E6930'] = 1.73205081
        v_['E6931'] = 3.00000000
        v_['E6932'] = 2.44948974
        v_['E6933'] = 2.64575131
        v_['E6934'] = 3.00000000
        v_['E6935'] = 3.46410162
        v_['E6936'] = 3.60555128
        v_['E6937'] = 2.64575131
        v_['E6938'] = 3.00000000
        v_['E6939'] = 3.60555128
        v_['E6940'] = 3.74165739
        v_['E6941'] = 3.00000000
        v_['E6942'] = 3.60555128
        v_['E6943'] = 4.47213595
        v_['E6944'] = 2.23606798
        v_['E6945'] = 1.73205081
        v_['E6946'] = 1.41421356
        v_['E6947'] = 2.00000000
        v_['E6948'] = 2.00000000
        v_['E6949'] = 2.44948974
        v_['E6950'] = 1.73205081
        v_['E6951'] = 2.23606798
        v_['E6952'] = 2.64575131
        v_['E6953'] = 1.73205081
        v_['E6954'] = 1.41421356
        v_['E6955'] = 1.00000000
        v_['E6956'] = 1.41421356
        v_['E6957'] = 1.73205081
        v_['E6958'] = 1.41421356
        v_['E6959'] = 1.41421356
        v_['E6960'] = 1.41421356
        v_['E6961'] = 2.00000000
        v_['E6962'] = 1.00000000
        v_['E6963'] = 1.00000000
        v_['E6964'] = 2.00000000
        v_['E6965'] = 1.73205081
        v_['E6966'] = 1.41421356
        v_['E6967'] = 1.00000000
        v_['E6968'] = 1.41421356
        v_['E6969'] = 1.00000000
        v_['E6970'] = 1.00000000
        v_['E6971'] = 1.00000000
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('C',ix_)
        pb.xnames=arrset(pb.xnames,iv,'C')
        [iv,ix_,_] = s2x_ii('L',ix_)
        pb.xnames=arrset(pb.xnames,iv,'L')
        [iv,ix_,_] = s2x_ii('A',ix_)
        pb.xnames=arrset(pb.xnames,iv,'A')
        [iv,ix_,_] = s2x_ii('B',ix_)
        pb.xnames=arrset(pb.xnames,iv,'B')
        [iv,ix_,_] = s2x_ii('I',ix_)
        pb.xnames=arrset(pb.xnames,iv,'I')
        [iv,ix_,_] = s2x_ii('S',ix_)
        pb.xnames=arrset(pb.xnames,iv,'S')
        [iv,ix_,_] = s2x_ii('X0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X0')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            v_['XOVERE'] = v_['EINV']*v_['X'+str(I)]
            [ig,ig_,_] = s2x_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
            iv = ix_['C']
            pbm.A[ig,iv] = v_['EINV']+pbm.A[ig,iv]
            iv = ix_['L']
            pbm.A[ig,iv] = v_['XOVERE']+pbm.A[ig,iv]
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
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            v_['YOVERE'] = v_['EINV']*v_['Y'+str(I)]
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],v_['YOVERE'])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('C' in ix_):
            pb.x0[ix_['C']] = 0.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['C']),0.0)
        if('L' in ix_):
            pb.x0[ix_['L']] = 0.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['L']),0.0)
        if('A' in ix_):
            pb.x0[ix_['A']] = 1.0
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['A']),1.0)
        if('B' in ix_):
            pb.x0[ix_['B']] = 0.05
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B']),0.05)
        if('I' in ix_):
            pb.x0[ix_['I']] = 597.076
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['I']),597.076)
        if('S' in ix_):
            pb.x0[ix_['S']] = 22.9096
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['S']),22.9096)
        if('X0' in ix_):
            pb.x0[ix_['X0']] = 24027.5
        else:
            pb.y0 = arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X0']),24027.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'B2BEXP', iet_)
        elftv = loaset(elftv,it,0,'A')
        elftv = loaset(elftv,it,1,'B')
        elftv = loaset(elftv,it,2,'I')
        elftv = loaset(elftv,it,3,'S')
        elftv = loaset(elftv,it,4,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            ename = 'B'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'B2BEXP')
            ielftype = arrset(ielftype, ie, iet_["B2BEXP"])
            vname = 'A'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='A')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='B')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = I
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='I')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'S'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='S')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X0'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['X'+str(I)])
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['MLOWER']),int(v_['MUPPER'])+1):
            v_['E'] = v_['E'+str(I)]
            v_['EINV'] = v_['ONE']/v_['E']
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,v_['EINV'])
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        pb.pbclass = "NOR2-MN-7-56"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar  = (
              arrset( pbm.efpar,0,np.sqrt(1.0e0/np.arctan(1.0e0)))    # this is TORPI)
        pbm.efpar = arrset( pbm.efpar,1,np.sqrt(0.5e0))    # this is ROOTP5
        return pbm

    @staticmethod
    def B2BEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        APB = EV_[0]+EV_[1]
        APB2 = APB*APB
        A2 = EV_[0]*EV_[0]
        B2 = EV_[1]*EV_[1]
        AB = EV_[0]*EV_[1]
        S2 = EV_[3]*EV_[3]
        S3 = EV_[3]*S2
        PI = 0.5e0*AB/APB
        P = EV_[2]*PI
        PAI = 0.5e0*EV_[1]/APB-0.5e0*AB/APB2
        PBI = 0.5e0*EV_[0]/APB-0.5e0*AB/APB2
        PA = PAI*EV_[2]
        loc_PB = PBI*EV_[2]
        PAB = EV_[2]*AB/APB**3
        PAA = -EV_[2]*EV_[1]/APB2+PAB
        PBB = -EV_[2]*EV_[0]/APB2+PAB
        XMY = pbm.elpar[iel_][0]-EV_[4]
        Z = XMY/EV_[3]
        ZY = -1.0e0/EV_[3]
        ZS = -XMY/EV_[3]**2
        ZSY = 1.0e0/EV_[3]**2
        ZSS = 2.0e0*XMY/EV_[3]**3
        R = np.exp(-0.5e0*Z**2)
        DR = -Z*R
        D2R = -R-Z*DR
        RS = DR*ZS
        RY = DR*ZY
        RSS = D2R*ZS*ZS+DR*ZSS
        RSY = D2R*ZS*ZY+DR*ZSY
        RYY = D2R*ZY*ZY
        AC = pbm.efpar[1]*(EV_[0]*EV_[3]+XMY/EV_[3])
        ACA = pbm.efpar[1]*EV_[3]
        ACS = pbm.efpar[1]*(EV_[0]-XMY/S2)
        ACY = -pbm.efpar[1]/EV_[3]
        ACAS = pbm.efpar[1]
        ACSS = 2.0e0*pbm.efpar[1]*XMY/S3
        ACSY = pbm.efpar[1]/S2
        BC = pbm.efpar[1]*(EV_[1]*EV_[3]+XMY/EV_[3])
        BCB = ACA
        BCS = pbm.efpar[1]*(EV_[1]-XMY/S2)
        BCY = ACY
        BCBS = pbm.efpar[1]
        BCSS = ACSS
        BCSY = ACSY
        QA = ERFC_EV_[3]CEV_[0]LED(AC)
        DQA = 2.0e0*AC*QA-pbm.efpar[0]
        D2QA = 2.0e0*(QA+AC*DQA)
        QAA = DQA*ACA
        QAS = DQA*ACS
        QAY = DQA*ACY
        QAAA = D2QA*ACA*ACA
        QAAS = D2QA*ACA*ACS+DQA*ACAS
        QAAY = D2QA*ACA*ACY
        QASS = D2QA*ACS*ACS+DQA*ACSS
        QASY = D2QA*ACS*ACY+DQA*ACSY
        QAYY = D2QA*ACY*ACY
        QB = ERFC_EV_[3]CEV_[0]LED(BC)
        DQB = 2.0e0*BC*QB-pbm.efpar[0]
        D2QB = 2.0e0*(QB+BC*DQB)
        QBB = DQB*BCB
        QBS = DQB*BCS
        QBY = DQB*BCY
        QBBB = D2QB*BCB*BCB
        QBBS = D2QB*BCB*BCS+DQB*BCBS
        QBBY = D2QB*BCB*BCY
        QBSS = D2QB*BCS*BCS+DQB*BCSS
        QBSY = D2QB*BCS*BCY+DQB*BCSY
        QBYY = D2QB*BCY*BCY
        T = QA+QB
        TA = QAA
        TB = QBB
        TS = QAS+QBS
        TY = QAY+QBY
        TAA = QAAA
        TAS = QAAS
        TAY = QAAY
        TBB = QBBB
        TBS = QBBS
        TBY = QBBY
        TSS = QASS+QBSS
        TSY = QASY+QBSY
        TYY = QAYY+QBYY
        f_   = P*T*R
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P*TA+PA*T)*R
            g_[1] = (P*TB+loc_PB*T)*R
            g_[2] = PI*T*R
            g_[3] = P*(T*RS+TS*R)
            g_[4] = P*(T*RY+TY*R)
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = R*(P*TAA+PAA*T+2.0e0*PA*TA)
                H_[0,1] = R*(PA*TB+loc_PB*TA+PAB*T)
                H_[1,0] = H_[0,1]
                H_[0,2] = (PI*TA+PAI*T)*R
                H_[2,0] = H_[0,2]
                H_[0,3] = (P*TA+PA*T)*RS+(P*TAS+PA*TS)*R
                H_[3,0] = H_[0,3]
                H_[0,4] = (P*TA+PA*T)*RY+(P*TAY+PA*TY)*R
                H_[4,0] = H_[0,4]
                H_[1,1] = R*(P*TBB+PBB*T+2.0e0*loc_PB*TB)
                H_[1,2] = (PI*TB+PBI*T)*R
                H_[2,1] = H_[1,2]
                H_[1,3] = (P*TB+loc_PB*T)*RS+(P*TBS+loc_PB*TS)*R
                H_[3,1] = H_[1,3]
                H_[1,4] = (P*TB+loc_PB*T)*RY+(P*TBY+loc_PB*TY)*R
                H_[4,1] = H_[1,4]
                H_[2,3] = PI*(T*RS+TS*R)
                H_[3,2] = H_[2,3]
                H_[2,4] = PI*(T*RY+TY*R)
                H_[4,2] = H_[2,4]
                H_[3,3] = P*(T*RSS+TSS*R+2.0e0*TS*RS)
                H_[3,4] = P*(T*RSY+TSY*R+TS*RY+TY*RS)
                H_[4,3] = H_[3,4]
                H_[4,4] = P*(T*RYY+TYY*R+2.0e0*TY*RY)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

