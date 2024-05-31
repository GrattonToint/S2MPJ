from s2xlib import *
class  LIN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LIN
#    *********
# 
#    A non-convex global optimization chemical equilibrium problem from the 
#    thesis of W.J. Lin.
#    It has a nonlinear objective and linear constraints.
# 
#    Source: illustrative example (section 4.6) in
#    C.M. McDonald and C.A. Floudas, "Global optimization for the phase 
#    and chemical equilibrium problem: application to the NRTL equation",
#    Computers & Chemical Engineering, (submitted), 1994.
# 
#    SIF input: Marcel Mongeau, 9 February 1994.
# 
#    classification = "OLR2-AY-4-2"
# 
#    PARAMETERS likely to be changed for different problems:
# 
#    Number of variable sets (# of phases)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LIN'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LIN'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['P'] = 2
        v_['C'] = 2
        v_['1'] = 1
        v_['2'] = 2
        v_['TAU'+str(int(v_['1']))+','+str(int(v_['2']))] = 3.00498
        v_['TAU'+str(int(v_['2']))+','+str(int(v_['1']))] = 4.69071
        v_['ALF'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.391965
        v_['ALF'+str(int(v_['2']))+','+str(int(v_['1']))] = 0.391965
        v_['INIT'+str(int(v_['1']))] = 0.5
        v_['INIT'+str(int(v_['2']))] = 0.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['C'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(K),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(K))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['C'])+1):
            for K in range(int(v_['1']),int(v_['P'])+1):
                [ig,ig_,_] = s2x_ii('MB'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'MB'+str(I))
                iv = ix_['X'+str(I)+','+str(K)]
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
        for I in range(int(v_['1']),int(v_['C'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['MB'+str(I)],float(v_['INIT'+str(I)]))
        v_['ZERO'] = 0.0
        v_['ONE'] = 1.0
        for I in range(int(v_['1']),int(v_['C'])+1):
            v_['ALF'+str(I)+','+str(I)] = v_['ZERO']
            v_['TAU'+str(I)+','+str(I)] = v_['ZERO']
        for I in range(int(v_['1']),int(v_['C'])+1):
            for J in range(int(v_['1']),int(v_['C'])+1):
                v_['MALF'] = -1.0*v_['ALF'+str(I)+','+str(J)]
                v_['PROD'] = v_['MALF']*v_['TAU'+str(I)+','+str(J)]
                v_['G'+str(I)+','+str(J)] = np.exp(v_['PROD'])
        for I in range(int(v_['1']),int(v_['C'])+1):
            v_['G'+str(I)+','+str(I)] = v_['ONE']
        for I in range(int(v_['1']),int(v_['C'])+1):
            for J in range(int(v_['1']),int(v_['C'])+1):
                v_['M'+str(I)+','+str(J)] = (v_['G'+str(I)+','+str(J)]*v_['TAU'+str(I)+
                     ','+str(J)])
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),1.e-12)
        pb.xupper = np.full((pb.n,1),+float('inf'))
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = v_['INIT1']
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['2']))]] = v_['INIT1']
        pb.xupper[ix_['X'+str(int(v_['2']))+','+str(int(v_['1']))]] = v_['INIT2']
        pb.xupper[ix_['X'+str(int(v_['2']))+','+str(int(v_['2']))]] = v_['INIT2']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = float(0.5)
        pb.x0[ix_['X'+str(int(v_['1']))+','+str(int(v_['2']))]] = float(0.0)
        pb.x0[ix_['X'+str(int(v_['2']))+','+str(int(v_['1']))]] = float(0.0)
        pb.x0[ix_['X'+str(int(v_['2']))+','+str(int(v_['2']))]] = float(0.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eXTAUG1', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftp = []
        elftp = loaset(elftp,it,0,'G11')
        elftp = loaset(elftp,it,1,'G12')
        elftp = loaset(elftp,it,2,'G21')
        elftp = loaset(elftp,it,3,'G22')
        elftp = loaset(elftp,it,4,'M11')
        elftp = loaset(elftp,it,5,'M12')
        elftp = loaset(elftp,it,6,'M21')
        elftp = loaset(elftp,it,7,'M22')
        [it,iet_,_] = s2x_ii( 'eXTAUG2', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftp = loaset(elftp,it,0,'G11')
        elftp = loaset(elftp,it,1,'G12')
        elftp = loaset(elftp,it,2,'G21')
        elftp = loaset(elftp,it,3,'G22')
        elftp = loaset(elftp,it,4,'M11')
        elftp = loaset(elftp,it,5,'M12')
        elftp = loaset(elftp,it,6,'M21')
        elftp = loaset(elftp,it,7,'M22')
        [it,iet_,_] = s2x_ii( 'eXLOGX', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'eXLOGXC', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y1')
        elftv = loaset(elftv,it,2,'Y2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for K in range(int(v_['1']),int(v_['P'])+1):
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXTAUG1')
            ielftype = arrset(ielftype, ie, iet_["eXTAUG1"])
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['1']))+','+str(K)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['2']))+','+str(K)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G11')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['1']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G12')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['1']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G21')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['2']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G22')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['2']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M11')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['1']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M12')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['1']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M21')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['2']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['1']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M22')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['2']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXTAUG2')
            ielftype = arrset(ielftype, ie, iet_["eXTAUG2"])
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['1']))+','+str(K)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'X'+str(int(v_['2']))+','+str(K)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G11')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['1']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G12')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['1']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G21')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['2']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G22')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['G'+str(int(v_['2']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M11')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['1']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M12')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['1']))+','+str(int(v_['2']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M21')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['2']))+','+str(int(v_['1']))])))
            ename = 'A'+str(int(v_['2']))+','+str(K)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='M22')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['M'+str(int(v_['2']))+','+str(int(v_['2']))])))
        for K in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['C'])+1):
                ename = 'B'+str(I)+','+str(K)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eXLOGX')
                ielftype = arrset(ielftype, ie, iet_["eXLOGX"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for K in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['C'])+1):
                ename = 'C'+str(I)+','+str(K)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eXLOGXC')
                ielftype = arrset(ielftype, ie, iet_["eXLOGXC"])
                vname = 'X'+str(I)+','+str(K)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['1']))+','+str(K)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['2']))+','+str(K)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,1.e-12,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['1']),int(v_['P'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['1']))+','+str(K)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['2']))+','+str(K)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
            for I in range(int(v_['1']),int(v_['C'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(K)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = "OLR2-AY-4-2"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXTAUG1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        INSUM1 = EV_[0]*pbm.elpar[iel_][0]+EV_[1]*pbm.elpar[iel_][2]
        INSUM2 = EV_[0]*pbm.elpar[iel_][1]+EV_[1]*pbm.elpar[iel_][3]
        RATIO1 = pbm.elpar[iel_][4]/INSUM1
        RATIO2 = pbm.elpar[iel_][5]/INSUM2
        TERM1 = EV_[0]*RATIO1
        TERM2 = EV_[1]*RATIO2
        SUM = TERM1+TERM2
        SQ1 = TERM1/INSUM1
        SQ2 = TERM2/INSUM2
        SQ11 = SQ1*pbm.elpar[iel_][0]
        SQ12 = SQ2*pbm.elpar[iel_][1]
        SQ21 = SQ1*pbm.elpar[iel_][2]
        SQ22 = SQ2*pbm.elpar[iel_][3]
        TRI1 = RATIO1-SQ11-SQ12
        TRI2 = RATIO2-SQ21-SQ22
        CUB1 = SQ11/INSUM1
        CUB2 = SQ12/INSUM2
        CUB11 = CUB1*pbm.elpar[iel_][0]
        CUB12 = CUB2*pbm.elpar[iel_][1]
        CUBM21 = CUB1*pbm.elpar[iel_][2]
        CUBM22 = CUB2*pbm.elpar[iel_][3]
        CUB21 = SQ21*pbm.elpar[iel_][2]/INSUM1
        CUB22 = SQ22*pbm.elpar[iel_][3]/INSUM2
        H1 = RATIO2-SQ22-2*SQ21
        H2 = pbm.elpar[iel_][5]*pbm.elpar[iel_][1]/INSUM2**2
        H3 = SQ1*pbm.elpar[iel_][2]**2/INSUM1
        H4 = pbm.elpar[iel_][5]*pbm.elpar[iel_][3]/INSUM2**2
        H5 = SQ2*pbm.elpar[iel_][3]**2/INSUM2
        f_   = EV_[0]*SUM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SUM+EV_[0]*TRI1
            g_[1] = EV_[0]*TRI2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2*(TRI1+EV_[0]*(-SQ11+CUB11+CUB12))
                H_[0,1] = H1+EV_[0]*(-H2+2*(CUBM21+CUBM22))
                H_[1,0] = H_[0,1]
                H_[1,1] = 2*EV_[0]*(H3-H4+H5)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXTAUG2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        INSUM1 = EV_[0]*pbm.elpar[iel_][0]+EV_[1]*pbm.elpar[iel_][2]
        INSUM2 = EV_[0]*pbm.elpar[iel_][1]+EV_[1]*pbm.elpar[iel_][3]
        RATIO1 = pbm.elpar[iel_][6]/INSUM1
        RATIO2 = pbm.elpar[iel_][7]/INSUM2
        TERM1 = EV_[0]*RATIO1
        TERM2 = EV_[1]*RATIO2
        SUM = TERM1+TERM2
        SQ1 = TERM1/INSUM1
        SQ2 = TERM2/INSUM2
        SQ11 = SQ1*pbm.elpar[iel_][0]
        SQ12 = SQ2*pbm.elpar[iel_][1]
        SQ21 = SQ1*pbm.elpar[iel_][2]
        SQ22 = SQ2*pbm.elpar[iel_][3]
        TRI1 = RATIO1-SQ11-SQ12
        TRI2 = RATIO2-SQ21-SQ22
        CUB1 = SQ11/INSUM1
        CUB2 = SQ12/INSUM2
        CUB11 = CUB1*pbm.elpar[iel_][0]
        CUB12 = CUB2*pbm.elpar[iel_][1]
        CUBM21 = CUB1*pbm.elpar[iel_][2]
        CUBM22 = CUB2*pbm.elpar[iel_][3]
        CUB21 = SQ21*pbm.elpar[iel_][2]/INSUM1
        CUB22 = SQ22*pbm.elpar[iel_][3]/INSUM2
        H1 = RATIO1-SQ11-2*SQ12
        H2 = pbm.elpar[iel_][6]*pbm.elpar[iel_][2]/INSUM1**2
        H3 = SQ1*pbm.elpar[iel_][0]**2/INSUM1
        H4 = pbm.elpar[iel_][6]*pbm.elpar[iel_][0]/INSUM1**2
        H5 = SQ2*pbm.elpar[iel_][1]**2/INSUM2
        f_   = EV_[1]*SUM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*TRI1
            g_[1] = SUM+EV_[1]*TRI2
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2*EV_[1]*(H3-H4+H5)
                H_[0,1] = H1+EV_[1]*(-H2+2*(CUBM21+CUBM22))
                H_[1,0] = H_[0,1]
                H_[1,1] = 2*(TRI2+EV_[1]*(-SQ22+CUB21+CUB22))
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LOGX = np.log(EV_[0])
        f_   = EV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX+1.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0/EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXLOGXC(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        U_[0,0] = U_[0,0]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        LOGX = np.log(IV_[1])
        f_   = IV_[0]*LOGX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LOGX
            g_[1] = IV_[0]/IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0/IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = -IV_[0]/IV_[1]**2
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

