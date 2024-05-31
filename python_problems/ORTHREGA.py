from s2xlib import *
class  ORTHREGA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGA
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) an ellipse to a set of points
#    in the plane.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "QQR2-AN-V-V"
# 
#    Number of levels in the generation of the data points
#    ( number of data points =     4**LEVELS
#      number of variables   = 2 * 4**LEVELS + 5
#      number of constraints =     4**LEVELS         )
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORTHREGA'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'ORTHREGA'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['LEVELS'] = int(2);  #  SIF file default value
        else:
            v_['LEVELS'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE LEVELS              3              $-PARAMETER n = 133    original value
# IE LEVELS              4              $-PARAMETER n = 517
# IE LEVELS              5              $-PARAMETER n = 2053
# IE LEVELS              6              $-PARAMETER n = 8197
# IE LEVELS              7              $-PARAMETER n = 32773
# IE LEVELS              8              $-PARAMETER n = 131077
        v_['A'] = 9.0
        v_['B'] = 6.0
        v_['CX'] = 0.5
        v_['CY'] = 0.5
        v_['1'] = 1
        v_['PI'] = 3.1415926535
        v_['-A'] = -1.0*v_['A']
        v_['-B'] = -1.0*v_['B']
        v_['NPTS'] = 1
        v_['XD'+str(int(v_['1']))] = v_['CX']
        v_['YD'+str(int(v_['1']))] = v_['CY']
        for I in range(int(v_['1']),int(v_['LEVELS'])+1):
            v_['NP'] = 0+v_['NPTS']
            for J in range(int(v_['1']),int(v_['NP'])+1):
                v_['XZ'+str(J)] = v_['XD'+str(J)]
                v_['YZ'+str(J)] = v_['YD'+str(J)]
            v_['NPTS'] = 0
            for J in range(int(v_['1']),int(v_['NP'])+1):
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['A']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['A']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['B']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['-B']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['-A']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['-A']
                v_['NPTS'] = 1+v_['NPTS']
                v_['XD'+str(int(v_['NPTS']))] = v_['XZ'+str(J)]+v_['-B']
                v_['YD'+str(int(v_['NPTS']))] = v_['YZ'+str(J)]+v_['B']
            v_['A'] = v_['A']/v_['PI']
            v_['B'] = v_['B']/v_['PI']
            v_['-A'] = v_['-A']/v_['PI']
            v_['-B'] = v_['-B']/v_['PI']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2x_ii('H11',ix_)
        pb.xnames=arrset(pb.xnames,iv,'H11')
        [iv,ix_,_] = s2x_ii('H12',ix_)
        pb.xnames=arrset(pb.xnames,iv,'H12')
        [iv,ix_,_] = s2x_ii('H22',ix_)
        pb.xnames=arrset(pb.xnames,iv,'H22')
        [iv,ix_,_] = s2x_ii('G1',ix_)
        pb.xnames=arrset(pb.xnames,iv,'G1')
        [iv,ix_,_] = s2x_ii('G2',ix_)
        pb.xnames=arrset(pb.xnames,iv,'G2')
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            [ig,ig_,_] = s2x_ii('OX'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('OY'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['Y'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('E'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'E'+str(I))
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['OX'+str(I)],float(v_['XD'+str(I)]))
            pbm.gconst = arrset(pbm.gconst,ig_['OY'+str(I)],float(v_['YD'+str(I)]))
            pbm.gconst = arrset(pbm.gconst,ig_['E'+str(I)],float(1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('H11' in ix_):
            pb.x0[ix_['H11']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['H11']),float(1.0)))
        if('H12' in ix_):
            pb.x0[ix_['H12']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['H12']),float(0.0)))
        if('H22' in ix_):
            pb.x0[ix_['H22']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['H22']),float(1.0)))
        if('G1' in ix_):
            pb.x0[ix_['G1']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['G1']),float(0.0)))
        if('G2' in ix_):
            pb.x0[ix_['G2']] = float(0.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['G2']),float(0.0)))
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            if('X'+str(I) in ix_):
                pb.x0[ix_['X'+str(I)]] = float(v_['XD'+str(I)])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)]),float(v_['XD'+str(I)])))
            if('Y'+str(I) in ix_):
                pb.x0[ix_['Y'+str(I)]] = float(v_['YD'+str(I)])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y'+str(I)]),float(v_['YD'+str(I)])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eHXX', iet_)
        elftv = loaset(elftv,it,0,'H')
        elftv = loaset(elftv,it,1,'X')
        [it,iet_,_] = s2x_ii( 'eHXY', iet_)
        elftv = loaset(elftv,it,0,'H')
        elftv = loaset(elftv,it,1,'X')
        elftv = loaset(elftv,it,2,'Y')
        [it,iet_,_] = s2x_ii( 'eGX', iet_)
        elftv = loaset(elftv,it,0,'G')
        elftv = loaset(elftv,it,1,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eHXX')
            ielftype = arrset(ielftype, ie, iet_["eHXX"])
            vname = 'H11'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='H')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eHXY')
            ielftype = arrset(ielftype, ie, iet_["eHXY"])
            vname = 'H12'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='H')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eHXX')
            ielftype = arrset(ielftype, ie, iet_["eHXX"])
            vname = 'H22'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='H')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'ED'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eGX')
            ielftype = arrset(ielftype, ie, iet_["eGX"])
            vname = 'G1'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='G')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'EE'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eGX')
            ielftype = arrset(ielftype, ie, iet_["eGX"])
            vname = 'G2'
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='G')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
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
        for I in range(int(v_['1']),int(v_['NPTS'])+1):
            ig = ig_['OX'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['OY'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            ig = ig_['E'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EB'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EC'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ED'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EE'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-2.0))
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
        pb.pbclass = "QQR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eHXX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[1]
            g_[1] = 2.0*EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EV_[1]+EV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = EV_[0]+EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eHXY(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eGX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
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
