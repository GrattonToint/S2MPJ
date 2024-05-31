from s2xlib import *
class  HS45(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS45
#    *********
# 
#    Source: problem 45 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph.L. Toint, October 1990.
# 
#    classification = "OBR2-AN-5-0"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS45'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS45'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 5
        v_['1'] = 1
        v_['-120'] = -120.0
        v_['PAR'] = 1.0/v_['-120']
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
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(-2.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            pb.xupper[ix_['X'+str(I)]] = v_['RI']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(2.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en5PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'E1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'en5PR')
        ielftype = arrset(ielftype, ie, iet_["en5PR"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,2.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,2.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,2.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,2.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,2.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['PAR']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-AN-5-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en5PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]*EV_[3]*EV_[4]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*EV_[3]*EV_[4]
            g_[1] = EV_[0]*EV_[2]*EV_[3]*EV_[4]
            g_[2] = EV_[0]*EV_[1]*EV_[3]*EV_[4]
            g_[3] = EV_[0]*EV_[1]*EV_[2]*EV_[4]
            g_[4] = EV_[0]*EV_[1]*EV_[2]*EV_[3]
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,1] = EV_[2]*EV_[3]*EV_[4]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*EV_[3]*EV_[4]
                H_[2,0] = H_[0,2]
                H_[0,3] = EV_[1]*EV_[2]*EV_[4]
                H_[3,0] = H_[0,3]
                H_[0,4] = EV_[1]*EV_[2]*EV_[3]
                H_[4,0] = H_[0,4]
                H_[1,2] = EV_[0]*EV_[3]*EV_[4]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*EV_[2]*EV_[4]
                H_[3,1] = H_[1,3]
                H_[1,4] = EV_[0]*EV_[2]*EV_[3]
                H_[4,1] = H_[1,4]
                H_[2,3] = EV_[0]*EV_[1]*EV_[4]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[0]*EV_[1]*EV_[3]
                H_[4,2] = H_[2,4]
                H_[3,4] = EV_[0]*EV_[1]*EV_[2]
                H_[4,3] = H_[3,4]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

