from s2xlib import *
class  CAMSHAPE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CAMSHAPE
#    *********
# 
#    Maximize the area of the valve opening for one rotation of a convex cam 
#    with constraints on the curvature and on the radius of the cam
# 
#    This is problem 4 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, November 2000
# 
#    classification = "LOR2-AN-V-V"
# 
#    The number of discretization points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CAMSHAPE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CAMSHAPE'
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
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# IE N                   400            $-PARAMETER
# IE N                   800            $-PARAMETER
        v_['RV'] = 1.0
        v_['RMAX'] = 2.0
        v_['RMIN'] = 1.0
        v_['RAV'] = v_['RMIN']+v_['RMAX']
        v_['RAV'] = 0.5*v_['RAV']
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['ALPHA'] = 1.5
        v_['N+1'] = 1+v_['N']
        v_['5(N+1)'] = 5*v_['N+1']
        v_['5(N+1)'] = float(v_['5(N+1)'])
        v_['DTHETA'] = 2.0*v_['PI']
        v_['DTHETA'] = v_['DTHETA']/v_['5(N+1)']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['PIRV'] = v_['PI']*v_['RV']
        v_['PIRV/N'] = v_['PIRV']/v_['RN']
        v_['-PIRV/N'] = -1.0*v_['PIRV/N']
        v_['CDTHETA'] = np.cos(v_['DTHETA'])
        v_['2CDTHETA'] = 2.0*v_['CDTHETA']
        v_['ADTHETA'] = v_['ALPHA']*v_['DTHETA']
        v_['-ADTHETA'] = -1.0*v_['ADTHETA']
        v_['2ADTHETA'] = 2.0*v_['ADTHETA']
        v_['-RMIN'] = -1.0*v_['RMIN']
        v_['-RMAX'] = -1.0*v_['RMAX']
        v_['-2RMAX'] = -2.0*v_['RMAX']
        v_['RMIN2'] = v_['RMIN']*v_['RMIN']
        v_['RMIN2CD'] = v_['RMIN']*v_['2CDTHETA']
        v_['RMAX2CD'] = v_['RMAX']*v_['2CDTHETA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('R'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'R'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('AREA',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['R'+str(I)]
            pbm.A[ig,iv] = float(v_['-PIRV/N'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('CO'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'CO'+str(I))
        [ig,ig_,_] = s2x_ii('E1',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E1')
        iv = ix_['R'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['-RMIN'])+pbm.A[ig,iv]
        iv = ix_['R'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(v_['RMIN2CD'])+pbm.A[ig,iv]
        v_['R'] = v_['RMIN2CD']-v_['RMIN']
        [ig,ig_,_] = s2x_ii('E2',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E2')
        iv = ix_['R'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(v_['R'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E3',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E3')
        iv = ix_['R'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['-RMAX'])+pbm.A[ig,iv]
        iv = ix_['R'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(v_['RMAX2CD'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('E4',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'E4')
        iv = ix_['R'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['-2RMAX'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CU'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CU'+str(int(v_['0'])))
        iv = ix_['R'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('CU'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CU'+str(I))
            iv = ix_['R'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['R'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('CU'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CU'+str(int(v_['N'])))
        iv = ix_['R'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['E2'],float(v_['RMIN2']))
        v_['R'] = v_['-ADTHETA']+v_['RMIN']
        pbm.gconst = arrset(pbm.gconst,ig_['CU'+str(int(v_['0']))],float(v_['R']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CU'+str(I)],float(v_['-ADTHETA']))
        v_['R'] = v_['-ADTHETA']-v_['RMAX']
        pbm.gconst = arrset(pbm.gconst,ig_['CU'+str(int(v_['N']))],float(v_['R']))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[legrps] = np.full((pb.nle,1),float('inf'))
        grange[gegrps] = np.full((pb.nge,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            grange = arrset(grange,ig_['CU'+str(I)],float(v_['2ADTHETA']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['R'+str(I)]] = v_['RMIN']
            pb.xupper[ix_['R'+str(I)]] = v_['RMAX']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.x0[ix_['R'+str(I)]] = float(v_['RAV'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ename = 'RA'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'R'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'R'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'RB'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'R'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'R'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'RA'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        ename = 'RA'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'R'+str(int(v_['N']))
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'RA'+str(int(v_['N']))
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'R'+str(int(v_['N-1']))
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'R2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
        ielftype = arrset(ielftype, ie, iet_["eSQR"])
        vname = 'R'+str(int(v_['N']))
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
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['CO'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RA'+str(int(v_['I+1']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['2CDTHETA']))
        ig = ig_['E1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RA'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['E3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RA'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        ig = ig_['E4']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['R2'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['2CDTHETA']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle)] = grange[legrps]
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        pb.cupper[np.arange(pb.nge)] = grange[gegrps]
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

