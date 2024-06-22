from s2mpjlib import *
class  TRAINF(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TRAINF
#    *********
# 
#    This is an optimal control problem.
#    The problem is to minimize the energy spent to move a train 
#    from the beginning of a flat track to its end in a given time.  The train
#    is slowed down by some drag (assumed to be quadratic in the the velocity).
#    The control variables are the acceleration force (UA) and the braking
#    force (UB) applied on the train.
# 
#    Source:
#    J. Kautsky and N. K. Nichols,
#    "OTEP-2: Optimal Train Energy Programme, mark 2",
#    Numerical Analysis Report NA/4/83,
#    Department of Mathematics, University of Reading, 1983.
# 
#    SIF input: N. Nichols and Ph. Toint, April 1993
# 
#    classification = "QQR2-MN-V-V"
# 
#    Problem variants
# 
#           Alternative values for the SIF file parameters:
# RE TIME                4.8            $-PARAMETER  travel time
# RE LENGTH              6.0            $-PARAMETER  length of track
# 
# RE TIME                2.0            $-PARAMETER  travel time
# RE LENGTH              2.0            $-PARAMETER  length of track
# 
# RE TIME                1.5            $-PARAMETER  travel time
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TRAINF'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['TIME'] = float(1.5);  #  SIF file default value
        else:
            v_['TIME'] = float(args[0])
# RE LENGTH              2.0            $-PARAMETER  length of track
        if nargin<2:
            v_['LENGTH'] = float(2);  #  SIF file default value
        else:
            v_['LENGTH'] = float(args[1])
# IE N                   11             $-PARAMETER
# IE N                   51             $-PARAMETER
# IE N                   101            $-PARAMETER     original value
# IE N                   201            $-PARAMETER
# IE N                   501            $-PARAMETER
# IE N                   1001           $-PARAMETER
        if nargin<3:
            v_['N'] = int(11);  #  SIF file default value
        else:
            v_['N'] = int(args[2])
# IE N                   5001           $-PARAMETER
# IE N                   10001          $-PARAMETER
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = v_['TIME']/v_['RN']
        v_['H/2'] = 0.5*v_['H']
        v_['-H'] = -1.0*v_['H']
        v_['-H/2'] = -1.0*v_['H/2']
        v_['UAMAX'] = 10.0
        v_['UBMIN'] = -2.0
        v_['VMAX'] = 10.0
        v_['A'] = 0.3
        v_['B'] = 0.14
        v_['C'] = 0.16
        v_['0'] = 0
        v_['1'] = 1
        v_['BH/2'] = v_['B']*v_['H/2']
        v_['1+BH/2'] = 1.0+v_['BH/2']
        v_['BH/2-1'] = -1.0+v_['BH/2']
        v_['-AH'] = v_['A']*v_['-H']
        v_['LENGTH/N'] = v_['LENGTH']/v_['RN']
        v_['CH/2'] = v_['C']*v_['H/2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('V'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('UA'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UA'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('UB'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UB'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('ENERGY',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('XEQ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'XEQ'+str(I))
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['V'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('VEQ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'VEQ'+str(I))
            iv = ix_['V'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['1+BH/2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)]
            pbm.A[ig,iv] = float(v_['BH/2-1'])+pbm.A[ig,iv]
            iv = ix_['UA'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            iv = ix_['UA'+str(I)]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            iv = ix_['UB'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
            iv = ix_['UB'+str(I)]
            pbm.A[ig,iv] = float(v_['-H/2'])+pbm.A[ig,iv]
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
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['VEQ'+str(I)],float(v_['-AH']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['V'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['V'+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['UA'+str(int(v_['0']))]] = v_['UAMAX']
        pb.xupper[ix_['UA'+str(int(v_['0']))]] = v_['UAMAX']
        pb.xlower[ix_['UB'+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['UB'+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pb.xlower[ix_['X'+str(I)]] = -float('Inf')
            pb.xupper[ix_['X'+str(I)]] = +float('Inf')
            pb.xlower[ix_['V'+str(I)]] = -float('Inf')
            pb.xupper[ix_['V'+str(I)]] = +float('Inf')
            pb.xlower[ix_['UA'+str(I)]] = 0.0
            pb.xupper[ix_['UA'+str(I)]] = v_['UAMAX']
            pb.xlower[ix_['UB'+str(I)]] = v_['UBMIN']
            pb.xupper[ix_['UB'+str(I)]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['N']))]] = v_['LENGTH']
        pb.xupper[ix_['X'+str(int(v_['N']))]] = v_['LENGTH']
        pb.xlower[ix_['V'+str(int(v_['N']))]] = 0.0
        pb.xupper[ix_['V'+str(int(v_['N']))]] = 0.0
        pb.xlower[ix_['UA'+str(int(v_['N']))]] = 0.0
        pb.xupper[ix_['UA'+str(int(v_['N']))]] = 0.0
        pb.xlower[ix_['UB'+str(int(v_['N']))]] = v_['UBMIN']
        pb.xupper[ix_['UB'+str(int(v_['N']))]] = v_['UBMIN']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['X'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['V'+str(int(v_['0']))]] = float(0.0)
        pb.x0[ix_['UA'+str(int(v_['0']))]] = float(v_['UAMAX'])
        pb.x0[ix_['UB'+str(int(v_['0']))]] = float(0.0)
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RI'] = float(I)
            v_['PI'] = v_['LENGTH/N']*v_['RI']
            pb.x0[ix_['X'+str(I)]] = float(v_['PI'])
            pb.x0[ix_['V'+str(I)]] = float(v_['LENGTH/N'])
            pb.x0[ix_['UA'+str(I)]] = float(0.0)
            pb.x0[ix_['UB'+str(I)]] = float(0.0)
        pb.x0[ix_['X'+str(int(v_['N']))]] = float(v_['LENGTH'])
        pb.x0[ix_['V'+str(int(v_['N']))]] = float(0.0)
        pb.x0[ix_['UA'+str(int(v_['N']))]] = float(0.0)
        pb.x0[ix_['UB'+str(int(v_['N']))]] = float(v_['UBMIN'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'UU')
        elftv = loaset(elftv,it,1,'VV')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'VVV')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            ename = 'VISQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'V'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='VVV')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'UV'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'UA'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='UU')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'V'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='VV')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['VEQ'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['VISQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CH/2']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['VISQ'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['CH/2']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['ENERGY']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['UV'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION            3.09751881012
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
        pb.pbclass = "QQR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

