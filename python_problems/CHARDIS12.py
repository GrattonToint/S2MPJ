from s2xlib import *
class  CHARDIS12(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHARDIS12
#    *********
# 
#    Distribution of charges on a round plate (2D)
# 
#    SIF input: R. Felkel, Jun 1999.
#               correction by S. Gratton & Ph. Toint, May 2024
#    modifield version of CHARDIS1 (formulation corrected)
# 
#    classification = "OQR2-AY-V-V"
# 
#    Number of positive (or negative) charges -> Number of variables 2*NP1
# 
#           Alternative values for the SIF file parameters:
# IE NP1                 5              $-PARAMETER
# IE NP1                 8              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHARDIS12'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CHARDIS12'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NP1'] = int(20);  #  SIF file default value
        else:
            v_['NP1'] = int(args[0])
# IE NP1                 50             $-PARAMETER
# IE NP1                 100            $-PARAMETER
# IE NP1                 200            $-PARAMETER
# IE NP1                 500            $-PARAMETER
# IE NP1                 1000           $-PARAMETER
# IE NP1                 2000           $-PARAMETER
# IE NP1                 5000           $-PARAMETER
        v_['R'] = 1.0
        v_['R2'] = v_['R']*v_['R']
        v_['N'] = -1+v_['NP1']
        v_['NReal'] = float(v_['N'])
        v_['NP1Real'] = float(v_['NP1'])
        v_['halfPI'] = np.arcsin(1.0)
        v_['PI'] = 2.0*v_['halfPI']
        v_['2PI'] = 4.0*v_['halfPI']
        v_['4PI'] = 8.0*v_['halfPI']
        v_['4PIqN'] = v_['4PI']/v_['NReal']
        v_['2PIqN'] = v_['2PI']/v_['NReal']
        v_['PIqN'] = v_['PI']/v_['NReal']
        v_['RqN'] = v_['R']/v_['NReal']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2x_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                [ig,ig_,_] = s2x_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            [ig,ig_,_] = s2x_ii('RES'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'RES'+str(I))
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
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['RES'+str(I)],float(v_['R2']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            pb.xlower[ix_['X'+str(I)]] = -float('Inf')
            pb.xupper[ix_['X'+str(I)]] = +float('Inf')
            pb.xlower[ix_['Y'+str(I)]] = -float('Inf')
            pb.xupper[ix_['Y'+str(I)]] = +float('Inf')
        pb.xlower[ix_['X1']] = v_['R']
        pb.xupper[ix_['X1']] = v_['R']
        pb.xlower[ix_['Y1']] = 0.0
        pb.xupper[ix_['Y1']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            v_['I-'] = -1+I
            v_['RealI-'] = float(v_['I-'])
            v_['RealNP1-I'] = v_['NP1Real']-v_['RealI-']
            v_['PHII-'] = v_['2PIqN']*v_['RealI-']
            v_['RI-'] = v_['RqN']*v_['RealNP1-I']
            v_['XST'] = np.cos(v_['PHII-'])
            v_['YST'] = np.sin(v_['PHII-'])
            v_['XS'] = v_['XST']*v_['RI-']
            v_['YS'] = v_['YST']*v_['RI-']
            pb.x0[ix_['X'+str(I)]] = float(v_['XS'])
            pb.x0[ix_['Y'+str(I)]] = float(v_['YS'])
        pb.x0[ix_['X1']] = float(v_['R'])
        pb.x0[ix_['Y1']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        [it,iet_,_] = s2x_ii( 'eDIFSQR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                ename = 'X'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eDIFSQR')
                ielftype = arrset(ielftype, ie, iet_["eDIFSQR"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Y'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eDIFSQR')
                ielftype = arrset(ielftype, ie, iet_["eDIFSQR"])
                vname = 'Y'+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            ename = 'RX'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'RY'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset(ielftype, ie, iet_["eSQR"])
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gREZIP',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NP1'])+1):
            v_['I+'] = 1+I
            for J in range(int(v_['I+']),int(v_['NP1'])+1):
                ig = ig_['O'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gREZIP')
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        for I in range(int(v_['2']),int(v_['NP1'])+1):
            ig = ig_['RES'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RX'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RY'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OQR2-AY-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(pbm,nargout,*args):

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

    @staticmethod
    def eDIFSQR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-EV_[1])*(EV_[0]-EV_[1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-EV_[1])
            g_[1] = -2.0*(EV_[0]-EV_[1])
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = 2.0
                H_[0,1] = -2.0
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gREZIP(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 1.0/GVAR_
        if nargout>1:
            g_ = -1.0/(GVAR_*GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0/(GVAR_*GVAR_*GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

