from s2mpjlib import *
class  HATFLDENE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HATFLDENE
#    *********
# 
#    An exponential fitting  test problem from the OPTIMA user manual.
# 
#    Source:
#    "The OPTIMA user manual (issue No.8, p. 37)",
#    Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.
# 
#    SIF input: Ph. Toint, May 1990.
#    Nonlinear-equations version of HATFLDE.SIF, Nick Gould, Jan 2020.
# 
#    classification = "NOR2-AN-3-21"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HATFLDENE'

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
        v_['NG'] = 21
        v_['T1'] = 0.3
        v_['T2'] = 0.35
        v_['T3'] = 0.4
        v_['T4'] = 0.45
        v_['T5'] = 0.5
        v_['T6'] = 0.55
        v_['T7'] = 0.6
        v_['T8'] = 0.65
        v_['T9'] = 0.7
        v_['T10'] = 0.75
        v_['T11'] = 0.8
        v_['T12'] = 0.85
        v_['T13'] = 0.9
        v_['T14'] = 0.95
        v_['T15'] = 1.0
        v_['T16'] = 1.05
        v_['T17'] = 1.1
        v_['T18'] = 1.15
        v_['T19'] = 1.2
        v_['T20'] = 1.25
        v_['T21'] = 1.3
        v_['Z1'] = 1.561
        v_['Z2'] = 1.473
        v_['Z3'] = 1.391
        v_['Z4'] = 1.313
        v_['Z5'] = 1.239
        v_['Z6'] = 1.169
        v_['Z7'] = 1.103
        v_['Z8'] = 1.04
        v_['Z9'] = 0.981
        v_['Z10'] = 0.925
        v_['Z11'] = 0.8721
        v_['Z12'] = 0.8221
        v_['Z13'] = 0.7748
        v_['Z14'] = 0.73
        v_['Z15'] = 0.6877
        v_['Z16'] = 0.6477
        v_['Z17'] = 0.6099
        v_['Z18'] = 0.5741
        v_['Z19'] = 0.5403
        v_['Z20'] = 0.5084
        v_['Z21'] = 0.4782
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2mpj_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2mpj_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NG'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
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
        for I in range(int(v_['1']),int(v_['NG'])+1):
            v_['-ZI'] = -1.0*v_['Z'+str(I)]
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['-ZI']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(1.0)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(-1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(-1.0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(0.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXPEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        [it,iet_,_] = s2mpj_ii( 'ePEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['NG'])+1):
            ename = 'A'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eXPEXP')
            ielftype = arrset(ielftype, ie, iet_["eXPEXP"])
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['T'+str(I)]))
            ename = 'B'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePEXP')
            ielftype = arrset(ielftype, ie, iet_["ePEXP"])
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
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
        for I in range(int(v_['1']),int(v_['NG'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               5.120377D-07
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-AN-3-21"
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
        EX = np.exp(pbm.elpar[iel_][0]*EV_[0])
        f_   = EX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*EX
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]*EX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXPEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EX = np.exp(pbm.elpar[iel_][0]*EV_[1])
        f_   = EV_[0]*EX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EX
            g_[1] = pbm.elpar[iel_][0]*EV_[0]*EX
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = pbm.elpar[iel_][0]*EX
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]*EV_[0]*EX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

