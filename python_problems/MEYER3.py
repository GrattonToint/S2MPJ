from s2xlib import *
class  MEYER3(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MEYER3
#    *********
#    A problem arising in the analysis of the resistance of a
#    thermistor, as formulated by Meyer.
# 
#    This function  is a nonlinear least squares with 16 groups.  Each
#    group has a nonlinear element.
# 
#    Source:  Problem 10 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley #29 (p. 73).
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-RN-3-0"
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MEYER3'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MEYER3'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['16'] = 16
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('X1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X1')
        xscale = arrset(xscale,iv,0.01)
        [iv,ix_,_] = s2x_ii('X2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X2')
        xscale = arrset(xscale,iv,1000.0)
        [iv,ix_,_] = s2x_ii('X3',ix_)
        pb.xnames=arrset(pb.xnames,iv,'X3')
        xscale = arrset(xscale,iv,100.0)
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['G1'],float(34780.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G2'],float(28610.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G3'],float(23650.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G4'],float(19630.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G5'],float(16370.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G6'],float(13720.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G7'],float(11540.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G8'],float(9744.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G9'],float(8261.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G10'],float(7030.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G11'],float(6005.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G12'],float(5147.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G13'],float(4427.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G14'],float(3820.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G15'],float(3307.0))
        pbm.gconst = arrset(pbm.gconst,ig_['G16'],float(2872.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X1']] = float(0.02)
        pb.x0[ix_['X2']] = float(4000.0)
        pb.x0[ix_['X3']] = float(250.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eGAUSS', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['16'])+1):
            v_['5I'] = 5*I
            v_['45+5I'] = 45+v_['5I']
            v_['TI'] = float(v_['45+5I'])
            ename = 'E'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eGAUSS')
            ielftype = arrset(ielftype, ie, iet_["eGAUSS"])
            vname = 'X1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='T')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TI']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            ig = ig_['G'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-RN-3-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eGAUSS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TPV3 = pbm.elpar[iel_][0]+EV_[2]
        EXPA = np.exp(EV_[1]/TPV3)
        V1EXPA = EV_[0]*EXPA
        TPV3SQ = TPV3*TPV3
        H22 = V1EXPA/TPV3SQ
        MG3 = -EV_[1]*H22
        HT = EV_[1]/TPV3SQ
        T33 = HT+2.0/TPV3
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
            g_[1] = V1EXPA/TPV3
            g_[2] = MG3
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EXPA/TPV3
                H_[1,0] = H_[0,1]
                H_[0,2] = -HT*EXPA
                H_[2,0] = H_[0,2]
                H_[1,1] = H22
                H_[1,2] = -H22+MG3/TPV3
                H_[2,1] = H_[1,2]
                H_[2,2] = -MG3*T33
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

