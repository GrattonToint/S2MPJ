from s2xlib import *
class  HIMMELBFNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBFNE
#    *********
# 
#    A 4 variables data fitting problems by Himmelblau.
# 
#    Source: problem 32 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    See Buckley#76 (p. 66)
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Nonlinear-equations version of HIMMELBF.SIF, Nick Gould, Jan 2020.
# 
#    classification = "NOR2-AN-4-7"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBFNE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HIMMELBFNE'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['A1'] = 0.0
        v_['A2'] = 0.000428
        v_['A3'] = 0.001000
        v_['A4'] = 0.001610
        v_['A5'] = 0.002090
        v_['A6'] = 0.003480
        v_['A7'] = 0.005250
        v_['B1'] = 7.391
        v_['B2'] = 11.18
        v_['B3'] = 16.44
        v_['B4'] = 16.20
        v_['B5'] = 22.20
        v_['B6'] = 24.02
        v_['B7'] = 31.32
        v_['1'] = 1
        v_['7'] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        [iv,ix_,_] = s2x_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        [iv,ix_,_] = s2x_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        [iv,ix_,_] = s2x_ii('X4',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X4')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            pbm.gscale = arrset(pbm.gscale,ig,float(0.0001))
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.full((ngrp,1),1.0)
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(2.7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(2.7)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(90.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X2']),float(90.0)))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(1500.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(1500.0)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(10.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X4']),float(10.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eHF', iet_)
        elftv = loaset(elftv,it,0,'XA')
        elftv = loaset(elftv,it,1,'XB')
        elftv = loaset(elftv,it,2,'XC')
        elftv = loaset(elftv,it,3,'XD')
        elftp = []
        elftp = loaset(elftp,it,0,'A')
        elftp = loaset(elftp,it,1,'B')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['7'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eHF')
                ielftype = arrset( ielftype,ie,iet_['eHF'])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XA')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XB')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XC')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XD')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(I)]))
            posep = find(elftp[ielftype[ie]],lambda x:x=='B')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['B'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-AN-4-7"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHF(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U  = (
              EV_[0]*EV_[0]+pbm.elpar[iel_][0]*EV_[1]*EV_[1]+pbm.elpar[iel_][0]*pbm.elpar[iel_][0]*EV_[2]*EV_[2])
        V = pbm.elpar[iel_][1]*(1.0+pbm.elpar[iel_][0]*EV_[3]*EV_[3])
        V2 = V*V
        AB = pbm.elpar[iel_][0]*pbm.elpar[iel_][1]
        A2 = pbm.elpar[iel_][0]*pbm.elpar[iel_][0]
        T = -4.0*AB/V2
        f_   = U/V
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*EV_[0]/V
            g_[1] = 2.0*pbm.elpar[iel_][0]*EV_[1]/V
            g_[2] = 2.0*A2*EV_[2]/V
            g_[3] = -2.0*AB*EV_[3]*U/V2
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 2.0/V
                H_[0,3] = T*EV_[3]*EV_[0]
                H_[3,0] = H_[0,3]
                H_[1,1] = 2.0*pbm.elpar[iel_][0]/V
                H_[1,3] = T*pbm.elpar[iel_][0]*EV_[3]*EV_[1]
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*A2/V
                H_[2,3] = T*pbm.elpar[iel_][0]*EV_[3]*EV_[2]
                H_[3,2] = H_[2,3]
                H_[3,3] = -2.0*AB*U/V2+8.0*(AB*EV_[3])**2*U/(V2*V)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

