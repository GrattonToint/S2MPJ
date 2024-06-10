from s2mpjlib import *
class  OSBORNE1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OSBORNE1
#    *********
# 
#    Osborne first problem in 5 variables. This is a nonlinear equation version
#    of problem OSBORNEA.
# 
#    This function  is a nonlinear least squares with 33 groups.  Each
#    group has 2 nonlinear elements and one linear element.
# 
#    Source:  Problem 17 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See alos Buckley#32 (p. 77).
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "NOR2-MN-5-33"
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OSBORNE1'

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
        v_['M'] = 33
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        [iv,ix_,_] = s2mpj_ii('X4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4')
        [iv,ix_,_] = s2mpj_ii('X5',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X5')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            iv = ix_['X1']
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
        pbm.gconst = arrset(pbm.gconst,ig_['G1'],float(0.844))
        pbm.gconst = arrset(pbm.gconst,ig_['G2'],float(0.908))
        pbm.gconst = arrset(pbm.gconst,ig_['G3'],float(0.932))
        pbm.gconst = arrset(pbm.gconst,ig_['G4'],float(0.936))
        pbm.gconst = arrset(pbm.gconst,ig_['G5'],float(0.925))
        pbm.gconst = arrset(pbm.gconst,ig_['G6'],float(0.908))
        pbm.gconst = arrset(pbm.gconst,ig_['G7'],float(0.881))
        pbm.gconst = arrset(pbm.gconst,ig_['G8'],float(0.850))
        pbm.gconst = arrset(pbm.gconst,ig_['G9'],float(0.818))
        pbm.gconst = arrset(pbm.gconst,ig_['G10'],float(0.784))
        pbm.gconst = arrset(pbm.gconst,ig_['G11'],float(0.751))
        pbm.gconst = arrset(pbm.gconst,ig_['G12'],float(0.718))
        pbm.gconst = arrset(pbm.gconst,ig_['G13'],float(0.685))
        pbm.gconst = arrset(pbm.gconst,ig_['G14'],float(0.658))
        pbm.gconst = arrset(pbm.gconst,ig_['G15'],float(0.628))
        pbm.gconst = arrset(pbm.gconst,ig_['G16'],float(0.603))
        pbm.gconst = arrset(pbm.gconst,ig_['G17'],float(0.580))
        pbm.gconst = arrset(pbm.gconst,ig_['G18'],float(0.558))
        pbm.gconst = arrset(pbm.gconst,ig_['G19'],float(0.538))
        pbm.gconst = arrset(pbm.gconst,ig_['G20'],float(0.522))
        pbm.gconst = arrset(pbm.gconst,ig_['G21'],float(0.506))
        pbm.gconst = arrset(pbm.gconst,ig_['G22'],float(0.490))
        pbm.gconst = arrset(pbm.gconst,ig_['G23'],float(0.478))
        pbm.gconst = arrset(pbm.gconst,ig_['G24'],float(0.467))
        pbm.gconst = arrset(pbm.gconst,ig_['G25'],float(0.457))
        pbm.gconst = arrset(pbm.gconst,ig_['G26'],float(0.448))
        pbm.gconst = arrset(pbm.gconst,ig_['G27'],float(0.438))
        pbm.gconst = arrset(pbm.gconst,ig_['G28'],float(0.431))
        pbm.gconst = arrset(pbm.gconst,ig_['G29'],float(0.424))
        pbm.gconst = arrset(pbm.gconst,ig_['G30'],float(0.420))
        pbm.gconst = arrset(pbm.gconst,ig_['G31'],float(0.414))
        pbm.gconst = arrset(pbm.gconst,ig_['G32'],float(0.411))
        pbm.gconst = arrset(pbm.gconst,ig_['G33'],float(0.406))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(0.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(0.5)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(1.5)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(1.5)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(-1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(-1.0)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(0.01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4']),float(0.01)))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(0.02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(0.02)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePEXP', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['I-1'] = -1+I
            v_['ITI'] = 10*v_['I-1']
            v_['MTI'] = float(v_['ITI'])
            v_['TI'] = -1.0*v_['MTI']
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePEXP')
            ielftype = arrset(ielftype, ie, iet_["ePEXP"])
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TI']))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePEXP')
            ielftype = arrset(ielftype, ie, iet_["ePEXP"])
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X5'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TI']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               5.46489D-05
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
        pb.pbclass = "NOR2-MN-5-33"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPA = np.exp(pbm.elpar[iel_][0]*EV_[1])
        V1EXPA = EV_[0]*EXPA
        f_   = V1EXPA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPA
            g_[1] = pbm.elpar[iel_][0]*V1EXPA
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]*EXPA
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]*V1EXPA
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

