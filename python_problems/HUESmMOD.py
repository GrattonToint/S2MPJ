from s2xlib import *
class  HUESmMOD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HUESmMOD
#    *********
# 
#    Source: An inverse problem from astronomy,
#    reformulated as a convex quadratic program by
#    S. P. Hestis, SIAM Review 34 (1992) pp. 642-647.
# 
#    SIF input: Nick Gould, January 1993.
#    improvements by: Ruediger Franke (Ruediger.Franke@RZ.TU-Ilmenau.DE)
# 
#    classification = "QLR2-MN-V-V"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HUESmMOD'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HUESmMOD'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['K'] = int(10);  #  SIF file default value
        else:
            v_['K'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE K                   100            $-PARAMETER
# IE K                   1000           $-PARAMETER    original value
# IE K                   5000           $-PARAMETER
# IE K                   10000          $-PARAMETER
        v_['1'] = 1
        v_['RANGE'] = 1.0
        v_['3.0'] = 3.0
        v_['5.0'] = 5.0
        v_['RK'] = float(v_['K'])
        v_['DELTAX'] = v_['RANGE']/v_['RK']
        v_['DELTAX2'] = v_['DELTAX']*v_['DELTAX']
        v_['DELTAX3'] = v_['DELTAX2']*v_['DELTAX']
        v_['DELTAX5'] = v_['DELTAX3']*v_['DELTAX2']
        v_['DELTAX3/3'] = v_['DELTAX3']/v_['3.0']
        v_['DELTAX5/5'] = v_['DELTAX5']/v_['5.0']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['K'])+1):
            [iv,ix_,_] = s2x_ii('M'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'M'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['K'])+1):
            [ig,ig_,_] = s2x_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI2'] = v_['RI']*v_['RI']
            v_['RI3'] = v_['RI2']*v_['RI']
            v_['RI5'] = v_['RI3']*v_['RI2']
            v_['RI-1'] = float(v_['I-1'])
            v_['RI-12'] = v_['RI-1']*v_['RI-1']
            v_['RI-13'] = v_['RI-12']*v_['RI-1']
            v_['RI-15'] = v_['RI-13']*v_['RI-12']
            v_['DIFF3'] = v_['RI3']-v_['RI-13']
            v_['DIFF5'] = v_['RI5']-v_['RI-15']
            v_['COEFF1'] = v_['DIFF3']*v_['DELTAX3/3']
            v_['COEFF2'] = v_['DIFF5']*v_['DELTAX5/5']
            [ig,ig_,_] = s2x_ii('E1',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E1')
            iv = ix_['M'+str(I)]
            pbm.A[ig,iv] = float(v_['COEFF1'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E2',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E2')
            iv = ix_['M'+str(I)]
            pbm.A[ig,iv] = float(v_['COEFF2'])+pbm.A[ig,iv]
        v_['RK'] = float(v_['K'])
        v_['1/RK'] = 1.0/v_['RK']
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
        pbm.gconst = arrset(pbm.gconst,ig_['E1'],float(1835.2))
        pbm.gconst = arrset(pbm.gconst,ig_['E2'],float(909.8))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        pb.y0 = np.full((pb.m,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U1')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['K'])+1):
            ename = 'E'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset( ielftype,ie,iet_['eSQ'])
            vname = 'M'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['K'])+1):
            ig = ig_['OBJ'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/RK']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "QLR2-MN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

