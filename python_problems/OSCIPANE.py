from s2xlib import *
class  OSCIPANE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OSCIPANE
#    *********
# 
#    An "oscillating path" problem due to Yuri Nesterov
# 
#    SIF input: Nick Gould, Dec 2006.
# 
#    classification = "NOR2-AN-V-V"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   5              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OSCIPANE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'OSCIPANE'
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
# IE N                   15             $-PARAMETER
# IE N                   25             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# RE RHO                 1.0            $-PARAMETER    Nesterov's original value
        if nargin<2:
            v_['RHO'] = float(500.0);  #  SIF file default value
        else:
            v_['RHO'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['ONE'] = 1.0
        v_['1/RHO'] = v_['ONE']/v_['RHO']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('Q1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'Q1')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('Q'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'Q'+str(I))
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/RHO']))
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
        pbm.gconst = arrset(pbm.gconst,ig_['Q1'],float(0.5))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['1']))]] = float(-1.0)
        for I in range(int(v_['2']),int(v_['N'])+1):
            pb.x0[ix_['X'+str(I)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eCHEB', iet_)
        elftv = loaset(elftv,it,0,'TAU')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ename = 'P'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eCHEB')
                ielftype = arrset( ielftype,ie,iet_['eCHEB'])
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='TAU')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['N'])+1):
            ig = ig_['Q'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
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
        pb.pbclass = "NOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCHEB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 2.0*EV_[0]**2-1.0
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 4.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

