from s2xlib import *
class  TRAINH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TRAINH
#    *********
# 
#    The problem is to minimize the energy spent to move a train 
#    from the beginning of a track to its end in a given time.  The train
#    is slowed down by some drag (assumed to be quadratic in the the velocity).
#    The track follows the slope of a hill.  The track geometry is given
#    by the equation
# 
#           1                    1  ns-1                          x - z_i
#    g(x) = - ( a_1 + a_{ns} ) + -- SUM ( s_{i+1} - s_i ) arctan( ------- )
#           2                    pi  1                              eps
# 
#    where the z_i are the breakpoints between sections of the track, and where 
#    the s_i are the "slopes" on these sections (eps is a regularization
#    parameter). Here we have a track of the overall shape
# 
#                      ______
#                     /      \      z0 = 0, z1 = 2, z2 = 4, z3 = 6
#                    /        \     s1 = 2, s2 = 0, s3 = -2
#                   /          \    eps = 0.05
# 
#    The control variables are the acceleration force (UA) and the braking
#    force (UB) applied on the train.
# 
#    Source: adapted from
#    J. Kautsky and N. K. Nichols,
#    "OTEP-2: Optimal Train Energy Programme, mark 2",
#    Numerical Analysis Report NA/4/83,
#    Department of Mathematics, University of Reading, 1983.
# 
#    SIF input: N. Nichols and Ph. Toint, April 1993
# 
#    classification = "QOR2-MN-V-V"
# 
#    Number of discretized points in the interval
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TRAINH'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'TRAINH'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(11);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   51             $-PARAMETER n=208, m=102
# IE N                   101            $-PARAMETER n=408, m=202  original value
# IE N                   201            $-PARAMETER n=808, m=402
# IE N                   501            $-PARAMETER n=2008, m=1002 
# IE N                   1001           $-PARAMETER n=4008, m=2002
# IE N                   5001           $-PARAMETER n=20008, m=10002
        if nargin<2:
            v_['TIME'] = float(4.8);  #  SIF file default value
        else:
            v_['TIME'] = float(args[1])
        if nargin<3:
            v_['LENGTH'] = float(6.0);  #  SIF file default value
        else:
            v_['LENGTH'] = float(args[2])
        if nargin<4:
            v_['NS'] = int(3);  #  SIF file default value
        else:
            v_['NS'] = int(args[3])
        if nargin<5:
            v_['Z1'] = float(2.0);  #  SIF file default value
        else:
            v_['Z1'] = float(args[4])
        if nargin<6:
            v_['Z2'] = float(4.0);  #  SIF file default value
        else:
            v_['Z2'] = float(args[5])
        if nargin<7:
            v_['S1'] = float(2.0);  #  SIF file default value
        else:
            v_['S1'] = float(args[6])
        if nargin<8:
            v_['S2'] = float(0.0);  #  SIF file default value
        else:
            v_['S2'] = float(args[7])
        if nargin<9:
            v_['S3'] = float(-2.0);  #  SIF file default value
        else:
            v_['S3'] = float(args[8])
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
        v_['EPS'] = 0.05
        v_['0'] = 0
        v_['1'] = 1
        v_['PI'] = 3.1415926535
        v_['NS-1'] = -1+v_['NS']
        v_['BH/2'] = v_['B']*v_['H/2']
        v_['1+BH/2'] = 1.0+v_['BH/2']
        v_['BH/2-1'] = -1.0+v_['BH/2']
        v_['-AH'] = v_['A']*v_['-H']
        v_['LENGTH/N'] = v_['LENGTH']/v_['RN']
        v_['CH/2'] = v_['C']*v_['H/2']
        v_['SUMS'] = v_['S'+str(int(v_['1']))]+v_['S'+str(int(v_['NS']))]
        v_['-AVS'] = -0.5*v_['SUMS']
        v_['-AVSH'] = v_['-AVS']*v_['H']
        v_['CNST'] = v_['-AH']+v_['-AVSH']
        v_['H/2PI'] = v_['H/2']/v_['PI']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('V'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('UA'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UA'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('UB'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'UB'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('ENERGY',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('XEQ'+str(I),ig_)
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
            [ig,ig_,_] = s2x_ii('VEQ'+str(I),ig_)
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
            pbm.gconst = arrset(pbm.gconst,ig_['VEQ'+str(I)],float(v_['CNST']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
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
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'UU')
        elftv = loaset(elftv,it,1,'VV')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'VVV')
        [it,iet_,_] = s2x_ii( 'eATAN', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftp = []
        elftp = loaset(elftp,it,0,'ZZ')
        elftp = loaset(elftp,it,1,'E')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            ename = 'VISQ'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'V'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='VVV')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for J in range(int(v_['1']),int(v_['NS-1'])+1):
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eATAN')
                ielftype = arrset(ielftype, ie, iet_["eATAN"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='ZZ')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['Z'+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='E')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['EPS']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'UV'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'UA'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='UU')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'V'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
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
            for J in range(int(v_['1']),int(v_['NS-1'])+1):
                v_['J+1'] = 1+J
                v_['DS'] = v_['S'+str(int(v_['J+1']))]-v_['S'+str(J)]
                v_['WJ'] = v_['DS']*v_['H/2PI']
                ig = ig_['VEQ'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WJ']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WJ']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['ENERGY']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['UV'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H']))
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
        pb.pbclass = "QOR2-MN-V-V"
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
    def eATAN(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DX = EV_[0]-pbm.elpar[iel_][0]
        E2 = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        f_   = np.arctan(DX/pbm.elpar[iel_][1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][1]/(E2+DX*DX)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -2.0*DX*pbm.elpar[iel_][1]/(E2+DX*DX)**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

