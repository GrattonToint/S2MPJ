from s2xlib import *
class  MISRA1B(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MISRA1B
#    *********
# 
#    NIST Data fitting problem MISRA1B given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 * (1-(1+b2*x/2)**(-2)) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Misra, D., NIST (1978).  
#      Dental Research Monomolecular Adsorption Study.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "NOR2-MN-2-14"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MISRA1B'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MISRA1B'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 14
        v_['N'] = 2
        v_['1'] = 1
        v_['X1'] = 77.6
        v_['X2'] = 114.9
        v_['X3'] = 141.1
        v_['X4'] = 190.8
        v_['X5'] = 239.9
        v_['X6'] = 289.0
        v_['X7'] = 332.8
        v_['X8'] = 378.4
        v_['X9'] = 434.8
        v_['X10'] = 477.3
        v_['X11'] = 536.8
        v_['X12'] = 593.1
        v_['X13'] = 689.1
        v_['X14'] = 760.0
        v_['Y1'] = 10.07
        v_['Y2'] = 14.73
        v_['Y3'] = 17.94
        v_['Y4'] = 23.93
        v_['Y5'] = 29.61
        v_['Y6'] = 35.18
        v_['Y7'] = 40.02
        v_['Y8'] = 44.82
        v_['Y9'] = 50.76
        v_['Y10'] = 55.05
        v_['Y11'] = 61.01
        v_['Y12'] = 66.40
        v_['Y13'] = 75.47
        v_['Y14'] = 81.78
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2x_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
            iv = ix_['B1']
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('B1' in ix_):
            pb.x0[ix_['B1']] = float(500.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B1']),float(500.0)))
        if('B2' in ix_):
            pb.x0[ix_['B2']] = float(0.0001)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B2']),float(0.0001)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE3', iet_)
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
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE3')
            ielftype = arrset(ielftype, ie, iet_["eE3"])
            vname = 'B1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
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
        pb.pbclass = "NOR2-MN-2-14"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R = 1.0+0.5*EV_[1]*pbm.elpar[iel_][0]
        R2 = R*R
        R3 = R*R2
        f_   = -EV_[0]/R2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -1.0/R2
            g_[1] = EV_[0]*pbm.elpar[iel_][0]/R3
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]/R3
                H_[1,0] = H_[0,1]
                H_[1,1] = -1.5*EV_[0]*pbm.elpar[iel_][0]**2/R2**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

