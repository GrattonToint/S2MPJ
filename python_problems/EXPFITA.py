from s2mpjlib import *
class  EXPFITA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EXPFITA
#    *********
# 
#    One sided rational approximation to the exponential function, as
#    described by Powell.
# 
#    Source:
#    M.J.D. Powell,
#    "A tolerant algorithm for linearly constrained optimization
#    calculations"'
#    Mathematical Programming 45(3), pp.561--562, 1989.
# 
#    SDIF input: Ph. Toint and N. Gould, May 1990.
# 
#    classification = "OLR2-AN-5-22"
# 
#    Number of fitting points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'EXPFITA'

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
        v_['R'] = 11
        v_['1'] = 1
        v_['5.0'] = 5.0
        v_['R-1'] = -1+v_['R']
        v_['RR-1'] = float(v_['R-1'])
        v_['5/R-1'] = v_['5.0']/v_['RR-1']
        for I in range(int(v_['1']),int(v_['R'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            v_['T'+str(I)] = v_['RI-1']*v_['5/R-1']
            v_['ET'+str(I)] = np.exp(v_['T'+str(I)])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('P0',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P0')
        [iv,ix_,_] = s2mpj_ii('P1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P1')
        [iv,ix_,_] = s2mpj_ii('P2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'P2')
        [iv,ix_,_] = s2mpj_ii('Q1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q1')
        [iv,ix_,_] = s2mpj_ii('Q2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'Q2')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['R'])+1):
            v_['TM5'] = -5.0+v_['T'+str(I)]
            v_['TM5SQ'] = v_['TM5']*v_['TM5']
            v_['QC1'] = v_['TM5']*v_['ET'+str(I)]
            v_['QC2'] = v_['TM5SQ']*v_['ET'+str(I)]
            v_['-QC1'] = -1.0*v_['QC1']
            v_['-QC2'] = -1.0*v_['QC2']
            v_['2T'] = v_['T'+str(I)]*v_['T'+str(I)]
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I))
            iv = ix_['P0']
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['P1']
            pbm.A[ig,iv] = float(v_['T'+str(I)])+pbm.A[ig,iv]
            iv = ix_['P2']
            pbm.A[ig,iv] = float(v_['2T'])+pbm.A[ig,iv]
            iv = ix_['Q1']
            pbm.A[ig,iv] = float(v_['-QC1'])+pbm.A[ig,iv]
            iv = ix_['Q2']
            pbm.A[ig,iv] = float(v_['-QC2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('B'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'B'+str(I))
            iv = ix_['Q1']
            pbm.A[ig,iv] = float(v_['TM5'])+pbm.A[ig,iv]
            iv = ix_['Q2']
            pbm.A[ig,iv] = float(v_['TM5SQ'])+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['R'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)],float(v_['ET'+str(I)]))
            pbm.gconst = arrset(pbm.gconst,ig_['B'+str(I)],float(-0.99999))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('P0' in ix_):
            pb.x0[ix_['P0']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P0']),float(1.0)))
        if('P1' in ix_):
            pb.x0[ix_['P1']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P1']),float(1.0)))
        if('P2' in ix_):
            pb.x0[ix_['P2']] = float(6.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['P2']),float(6.0)))
        if('Q1' in ix_):
            pb.x0[ix_['Q1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q1']),float(0.0)))
        if('Q2' in ix_):
            pb.x0[ix_['Q2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q2']),float(0.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eFIT', iet_)
        elftv = loaset(elftv,it,0,'P0')
        elftv = loaset(elftv,it,1,'P1')
        elftv = loaset(elftv,it,2,'P2')
        elftv = loaset(elftv,it,3,'Q1')
        elftv = loaset(elftv,it,4,'Q2')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['R'])+1):
            ename = 'F'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eFIT')
            ielftype = arrset(ielftype, ie, iet_["eFIT"])
            vname = 'P0'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P0')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Q1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Q1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Q2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Q2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['T'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['R'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
        pass
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-AN-5-22"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eFIT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TM5 = pbm.elpar[iel_][0]-5.0
        TM5SQ = TM5*TM5
        T2 = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]
        ET = np.exp(pbm.elpar[iel_][0])
        QT = 1.0+EV_[3]*TM5+EV_[4]*TM5SQ
        ETQT = ET*QT
        ETQT2 = ETQT*QT
        ETQT3 = ETQT2*QT
        PT = EV_[0]+EV_[1]*pbm.elpar[iel_][0]+EV_[2]*T2
        F = PT/ETQT-1.0
        TWOF = F+F
        DFDP0 = 1.0/ETQT
        DFDP1 = pbm.elpar[iel_][0]/ETQT
        DFDP2 = T2/ETQT
        DFDQ1 = -PT*TM5/ETQT2
        DFDQ2 = -PT*TM5SQ/ETQT2
        D2P0Q1 = -TM5/ETQT2
        D2P0Q2 = -TM5SQ/ETQT2
        D2P1Q1 = -pbm.elpar[iel_][0]*TM5/ETQT2
        D2P1Q2 = -pbm.elpar[iel_][0]*TM5SQ/ETQT2
        D2P2Q1 = -T2*TM5/ETQT2
        D2P2Q2 = -T2*TM5SQ/ETQT2
        D2Q1Q1 = 2.0*PT*TM5SQ/ETQT3
        D2Q1Q2 = 2.0*PT*TM5SQ*TM5/ETQT3
        D2Q2Q2 = 2.0*PT*TM5SQ*TM5SQ/ETQT3
        f_   = F*F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TWOF*DFDP0
            g_[1] = TWOF*DFDP1
            g_[2] = TWOF*DFDP2
            g_[3] = TWOF*DFDQ1
            g_[4] = TWOF*DFDQ2
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = 2.0*DFDP0*DFDP0
                H_[0,1] = 2.0*DFDP0*DFDP1
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0*DFDP0*DFDP2
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*DFDP1*DFDP1
                H_[1,2] = 2.0*DFDP1*DFDP2
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*DFDP2*DFDP2
                H_[0,3] = TWOF*D2P0Q1+2.0*DFDP0*DFDQ1
                H_[3,0] = H_[0,3]
                H_[0,4] = TWOF*D2P0Q2+2.0*DFDP0*DFDQ2
                H_[4,0] = H_[0,4]
                H_[1,3] = TWOF*D2P1Q1+2.0*DFDP1*DFDQ1
                H_[3,1] = H_[1,3]
                H_[1,4] = TWOF*D2P1Q2+2.0*DFDP1*DFDQ2
                H_[4,1] = H_[1,4]
                H_[2,3] = TWOF*D2P2Q1+2.0*DFDP2*DFDQ1
                H_[3,2] = H_[2,3]
                H_[2,4] = TWOF*D2P2Q2+2.0*DFDP2*DFDQ2
                H_[4,2] = H_[2,4]
                H_[3,3] = TWOF*D2Q1Q1+2.0*DFDQ1*DFDQ1
                H_[3,4] = TWOF*D2Q1Q2+2.0*DFDQ1*DFDQ2
                H_[4,3] = H_[3,4]
                H_[4,4] = TWOF*D2Q2Q2+2.0*DFDQ2*DFDQ2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

