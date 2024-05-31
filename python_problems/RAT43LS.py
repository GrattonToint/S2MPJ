from s2xlib import *
class  RAT43LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RAT43LS
#    *********
# 
#    NIST Data fitting problem RAT43.
# 
#    Fit: y = b1 / ((1+exp[b2-b3*x])**(1/b4)) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Ratkowsky, D.A. (1983).  
#      Nonlinear Regression Modeling.
#      New York, NY:  Marcel Dekker, pp. 62 and 88.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "SUR2-MN-4-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RAT43LS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'RAT43LS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 15
        v_['N'] = 4
        v_['1'] = 1
        v_['X1'] = 1.0
        v_['X2'] = 2.0
        v_['X3'] = 3.0
        v_['X4'] = 4.0
        v_['X5'] = 5.0
        v_['X6'] = 6.0
        v_['X7'] = 7.0
        v_['X8'] = 8.0
        v_['X9'] = 9.0
        v_['X10'] = 10.0
        v_['X11'] = 11.0
        v_['X12'] = 12.0
        v_['X13'] = 13.0
        v_['X14'] = 14.0
        v_['X15'] = 15.0
        v_['Y1'] = 16.08
        v_['Y2'] = 33.83
        v_['Y3'] = 65.80
        v_['Y4'] = 97.20
        v_['Y5'] = 191.55
        v_['Y6'] = 326.20
        v_['Y7'] = 386.87
        v_['Y8'] = 520.53
        v_['Y9'] = 590.03
        v_['Y10'] = 651.92
        v_['Y11'] = 724.93
        v_['Y12'] = 699.56
        v_['Y13'] = 689.96
        v_['Y14'] = 637.56
        v_['Y15'] = 717.41
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
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
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
        pb.x0[ix_['B1']] = float(100.0)
        pb.x0[ix_['B2']] = float(10.0)
        pb.x0[ix_['B3']] = float(1.0)
        pb.x0[ix_['B4']] = float(1.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
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
            pbm.elftype = arrset(pbm.elftype,ie,'eE')
            ielftype = arrset(ielftype, ie, iet_["eE"])
            vname = 'B1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['X'+str(I)]))
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
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-MN-4-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V2MV3X = EV_[1]-EV_[2]*pbm.elpar[iel_][0]
        V4INV = 1.0/EV_[3]
        V4INVP = V4INV+1.0
        E = np.exp(V2MV3X)
        E2 = E*E
        EP1 = E+1.0
        EP1L = np.log(EP1)
        EP14 = EP1**V4INV
        EP14P1 = EP1**V4INVP
        EP14P2 = EP1**(V4INV+2.0)
        VE = EV_[3]*EP14P1
        VE2 = EV_[3]*EP14P2
        V42EPP = EP14*EV_[3]**2
        V42EP2 = EP14P1*EV_[3]**2
        V42EP3 = EP14P1*EV_[3]**3
        f_   = EV_[0]/EP14
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/EP14
            g_[1] = -EV_[0]*E/VE
            g_[2] = EV_[0]*pbm.elpar[iel_][0]*E/VE
            g_[3] = EV_[0]*EP1L/V42EPP
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = -E/VE
                H_[1,0] = H_[0,1]
                H_[0,2] = pbm.elpar[iel_][0]*E/VE
                H_[2,0] = H_[0,2]
                H_[0,3] = EP1L/V42EPP
                H_[3,0] = H_[0,3]
                H_[1,1] = EV_[0]*(E2*V4INVP/VE2-E/VE)
                H_[1,2] = EV_[0]*pbm.elpar[iel_][0]*(E/VE-E2*V4INVP/VE2)
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*E*(1.0/V42EP2-EP1L/V42EP3)
                H_[3,1] = H_[1,3]
                H_[2,2] = EV_[0]*pbm.elpar[iel_][0]**2*(E2*V4INVP/VE2-E/VE)
                H_[2,3] = EV_[0]*pbm.elpar[iel_][0]*E*(EP1L/V42EP3-1.0/V42EP2)
                H_[3,2] = H_[2,3]
                H_[3,3] = (EV_[0]/EP14)*(EP1L**2/EV_[3]**4-2.0*EP1L/EV_[3]**3)
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

