from s2xlib import *
class  COOLHANS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : COOLHANS
#    *********
# 
#    A problem arising from the analysis of a Cooley-Hansen economy with
#    loglinear approximation.  The problem is to solve the matrix equation
#                  A * X * X + B * X + C = 0
#    where A, B and C are known N times N matrices and X an unknown matrix
#    of matching dimension.  The instance considered here has N = 3.
# 
#    Source:
#    S. Ceria, private communication, 1995.
# 
#    SIF input: Ph. Toint, Feb 1995.
# 
#    classification = "NQR2-RN-9-9"
# 
#    order of the matrix equation
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'COOLHANS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'COOLHANS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 3
        v_['A1,1'] = 0.0
        v_['A2,1'] = 0.13725e-6
        v_['A3,1'] = 0.0
        v_['A1,2'] = 0.0
        v_['A2,2'] = 937.62
        v_['A3,2'] = 0.0
        v_['A1,3'] = 0.0
        v_['A2,3'] = -42.207
        v_['A3,3'] = 0.0
        v_['B1,1'] = 0.0060893
        v_['B2,1'] = 0.13880e-6
        v_['B3,1'] = -0.13877e-6
        v_['B1,2'] = -44.292
        v_['B2,2'] = -1886.0
        v_['B3,2'] = 42.362
        v_['B1,3'] = 2.0011
        v_['B2,3'] = 42.362
        v_['B3,3'] = -2.0705
        v_['C1,1'] = 0.0
        v_['C2,1'] = 0.0
        v_['C3,1'] = 0.0
        v_['C1,2'] = 44.792
        v_['C2,2'] = 948.21
        v_['C3,2'] = -42.684
        v_['C1,3'] = 0.0
        v_['C2,3'] = 0.0
        v_['C3,3'] = 0.0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    [ig,ig_,_] = s2x_ii('G'+str(K)+','+str(L),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'G'+str(K)+','+str(L))
                    iv = ix_['X'+str(M)+','+str(L)]
                    pbm.A[ig,iv] = float(v_['B'+str(K)+','+str(M)])+pbm.A[ig,iv]
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
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                v_['-C'] = -1.0*v_['C'+str(K)+','+str(L)]
                pbm.gconst = arrset(pbm.gconst,ig_['G'+str(K)+','+str(L)],float(v_['-C']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for K in range(int(v_['1']),int(v_['N'])+1):
            for L in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'E'+str(K)+','+str(M)+','+str(L)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                    ielftype = arrset(ielftype, ie, iet_["en2PR"])
                    pb.x0 = np.zeros((pb.n,1))
                    vname = 'X'+str(K)+','+str(M)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(M)+','+str(L)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['1']))+','+str(L)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['1']))+','+str(P)])))
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['2']))+','+str(L)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['2']))+','+str(P)])))
        for L in range(int(v_['1']),int(v_['N'])+1):
            for P in range(int(v_['1']),int(v_['N'])+1):
                for M in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['G'+str(int(v_['3']))+','+str(L)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['E'+str(P)+','+str(M)+','+str(L)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['A'+str(int(v_['3']))+','+str(P)])))
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
        pb.pbclass = "NQR2-RN-9-9"
        pb.x0          = np.zeros((pb.n,1))
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

