from s2xlib import *
class  HAGER4(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HAGER4
#    *********
# 
#    A nonlinear optimal control problem, by W. Hager.
# 
#    NOTE: The solution for x given in the article below by Hager has
#    a typo. On the interval [1/2, 1], x(t) = (exp(2t) + exp(t))/d. In
#    other words, the minus sign in the article should be a plus sign.
# 
#    Source: problem P4 in
#    W.W. Hager,
#    "Multiplier Methods for Nonlinear Optimal Control",
#    SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "OLR2-AN-V-V"
# 
#    Number of discretized points in [0,1]
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HAGER4'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HAGER4'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   2500           $-PARAMETER
# IE N                   5000           $-PARAMETER
        v_['1/H'] = float(v_['N'])
        v_['H'] = 1.0/v_['1/H']
        v_['H/2'] = 0.5*v_['H']
        v_['1/H-1'] = -1.0+v_['1/H']
        v_['-1/H'] = -1.0*v_['1/H']
        v_['1/HSQ'] = v_['1/H']*v_['1/H']
        v_['1/2HSQ'] = 0.5*v_['1/HSQ']
        v_['0'] = 0
        v_['1'] = 1
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'+str(I)] = v_['RI']*v_['H']
            v_['-2TI'] = -2.0*v_['T'+str(I)]
            v_['Z'+str(I)] = np.exp(v_['-2TI'])
        for I in range(int(v_['0']),int(v_['1'])+1):
            v_['A'+str(I)] = -0.5*v_['Z'+str(I)]
            v_['TI+1/2'] = 0.5+v_['T'+str(I)]
            v_['B'+str(I)] = v_['A'+str(I)]*v_['TI+1/2']
            v_['TISQ'] = v_['T'+str(I)]*v_['T'+str(I)]
            v_['TIETC'] = v_['TISQ']+v_['TI+1/2']
            v_['C'+str(I)] = v_['A'+str(I)]*v_['TIETC']
        v_['DA'] = v_['A'+str(int(v_['1']))]-v_['A'+str(int(v_['0']))]
        v_['SCDA'] = 0.5*v_['DA']
        v_['DB'] = v_['B'+str(int(v_['1']))]-v_['B'+str(int(v_['0']))]
        v_['SCDB'] = v_['DB']*v_['1/H']
        v_['DC'] = v_['C'+str(int(v_['1']))]-v_['C'+str(int(v_['0']))]
        v_['SCDC'] = v_['DC']*v_['1/2HSQ']
        v_['E'] = np.exp(1.0)
        v_['3E'] = 3.0*v_['E']
        v_['1+3E'] = 1.0+v_['3E']
        v_['1-E'] = 1.0-v_['E']
        v_['2-2E'] = 2.0*v_['1-E']
        v_['XX0'] = v_['1+3E']/v_['2-2E']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2x_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'S'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['1/H-1'])+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            v_['ETI'] = np.exp(v_['T'+str(I)])
            v_['-ETI'] = -1.0*v_['ETI']
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['-ETI'])+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['X'+str(int(v_['0']))]] = v_['XX0']
        pb.xupper[ix_['X'+str(int(v_['0']))]] = v_['XX0']
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['0']))]] = float(v_['XX0'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eELT', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'D')
        elftp = loaset(elftp,it,1,'E')
        elftp = loaset(elftp,it,2,'F')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ename = 'EL'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eELT')
            ielftype = arrset(ielftype, ie, iet_["eELT"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['DD'] = v_['SCDA']*v_['Z'+str(int(v_['I-1']))]
            v_['EE'] = v_['SCDB']*v_['Z'+str(int(v_['I-1']))]
            v_['FF'] = v_['SCDC']*v_['Z'+str(int(v_['I-1']))]
            posep = find(elftp[ielftype[ie]],lambda x:x=='D')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['DD']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='E')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['EE']))
            posep = find(elftp[ielftype[ie]],lambda x:x=='F')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['FF']))
            ename = 'U'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EL'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['U'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H/2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
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
        pb.pbclass = "OLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eELT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_    = (
              pbm.elpar[iel_][0]*EV_[0]*EV_[0]+pbm.elpar[iel_][1]*EV_[0]*(EV_[1]-EV_[0])+pbm.elpar[iel_][2]*(EV_[1]-EV_[0])**2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (2.0*pbm.elpar[iel_][0]*EV_[0]+pbm.elpar[iel_][1]*(EV_[1]-2.0*EV_[0])-
                 2.0*pbm.elpar[iel_][2]*(EV_[1]-EV_[0]))
            g_[1] = pbm.elpar[iel_][1]*EV_[0]+2.0*pbm.elpar[iel_][2]*(EV_[1]-EV_[0])
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0*(pbm.elpar[iel_][0]-pbm.elpar[iel_][1]+pbm.elpar[iel_][2])
                H_[0,1] = pbm.elpar[iel_][1]-2.0*pbm.elpar[iel_][2]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*pbm.elpar[iel_][2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

