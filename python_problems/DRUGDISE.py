from s2mpjlib import *
class  DRUGDISE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRUGDISE
#    *********
# 
#    This is a variant of the drug displacement problem DRUGDIS where the
#    state equations have been Expanded in term of more intermediate
#    functions, each one of them being less nonlinear.
# 
#    The problem is based on the kinetic model of Aarons and Rowland which
#    simulates the interaction of the two drugs (warfarin and phenylnutazone)
#    in a patient bloodstream.  The state variable are the concentrations of
#    unbound warfarin (w) and phenylbutazone (p).  The problem is to control
#    the rate of injection (u) of the pain-killing phenylbutazone so that both
#    drugs reach a specified steady-state in minimum time and the concentration
#    of warfarin does not rise above a toxicity level.
# 
#    The problem is discretized using the trapeziodal rule.  It is non-convex.
# 
#    Source:
#    H. Maurer and M. Wiegand,
#    "Numerical solution of a drug displacement problem with bounded state
#    variables",
#    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "LOR2-MY-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#           Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=63, m=50 
# IE NI                  100            $-PARAMETER n=603, m=500   original value
# IE NI                  100            $-PARAMETER n=6003, m=5000 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DRUGDISE'

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
            v_['NI'] = int(10);  #  SIF file default value
        else:
            v_['NI'] = int(args[0])
        if nargin<2:
            v_['TOXIC'] = float(0.026);  #  SIF file default value
        else:
            v_['TOXIC'] = float(args[1])
        if nargin<3:
            v_['WSS'] = float(0.02);  #  SIF file default value
        else:
            v_['WSS'] = float(args[2])
        if nargin<4:
            v_['UMAX'] = float(8.0);  #  SIF file default value
        else:
            v_['UMAX'] = float(args[3])
        if nargin<5:
            v_['PSTART'] = float(0.0);  #  SIF file default value
        else:
            v_['PSTART'] = float(args[4])
        if nargin<6:
            v_['PFINAL'] = float(2.0);  #  SIF file default value
        else:
            v_['PFINAL'] = float(args[5])
        if nargin<7:
            v_['Z'] = float(46.4);  #  SIF file default value
        else:
            v_['Z'] = float(args[6])
        v_['AVP'] = v_['PSTART']+v_['PFINAL']
        v_['AVP'] = 0.5*v_['AVP']
        v_['-Z'] = -1.0*v_['Z']
        v_['-ZZ'] = v_['Z']*v_['-Z']
        v_['NI-1'] = -1+v_['NI']
        v_['RNI'] = float(v_['NI'])
        v_['-1/NI'] = -1.0/v_['RNI']
        v_['-Z/NI'] = v_['Z']*v_['-1/NI']
        v_['0'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('TF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TF')
        xscale = arrset(xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I))
            xscale = arrset(xscale,iv,0.02)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('A'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'A'+str(I))
            xscale = arrset(xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
            xscale = arrset(xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'C'+str(I))
            xscale = arrset(xscale,iv,0.0000001)
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('TFINAL',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['TF']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(100.0))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('EW'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EW'+str(I))
            iv = ix_['W'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(0.02))
            [ig,ig_,_] = s2mpj_ii('EP'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EP'+str(I))
            iv = ix_['P'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['P'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('EA'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EA'+str(I))
            iv = ix_['A'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['P'+str(I)]
            pbm.A[ig,iv] = float(v_['-Z'])+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(200.0))
            [ig,ig_,_] = s2mpj_ii('EB'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EB'+str(I))
            iv = ix_['B'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['-Z'])+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(200.0))
            [ig,ig_,_] = s2mpj_ii('EC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EC'+str(I))
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
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['EA'+str(I)],float(232.0))
            pbm.gconst = arrset(pbm.gconst,ig_['EB'+str(I)],float(232.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            pb.xlower[ix_['C'+str(I)]] = -float('Inf')
            pb.xupper[ix_['C'+str(I)]] = +float('Inf')
        pb.xlower[ix_['TF']] = 200.0
        for I in range(int(v_['0']),int(v_['NI'])+1):
            pb.xupper[ix_['W'+str(I)]] = v_['TOXIC']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            pb.xupper[ix_['U'+str(I)]] = v_['UMAX']
        pb.xlower[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        pb.xupper[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        pb.xlower[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        pb.xupper[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        pb.xlower[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        pb.xupper[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        pb.xlower[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        pb.xupper[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        v_['2W/10'] = 0.2*v_['WSS']
        v_['2P/10'] = 0.2*v_['AVP']
        v_['2(W+P)/10'] = v_['2W/10']+v_['2P/10']
        v_['D'] = 1.0+v_['2(W+P)/10']
        v_['DD'] = v_['D']*v_['D']
        v_['ZP'] = v_['AVP']*v_['Z']
        v_['ZW'] = v_['WSS']*v_['Z']
        v_['AA'] = v_['DD']+v_['ZP']
        v_['AA'] = 232.0+v_['AA']
        v_['BB'] = v_['DD']+v_['ZW']
        v_['BB'] = 232.0+v_['BB']
        v_['AB'] = v_['AA']*v_['BB']
        v_['WP'] = v_['WSS']*v_['AVP']
        v_['-ZZWP'] = v_['WP']*v_['-ZZ']
        v_['CD'] = v_['AB']+v_['-ZZWP']
        v_['CC'] = v_['DD']/v_['CD']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            pb.x0[ix_['W'+str(I)]] = float(v_['WSS'])
            pb.x0[ix_['P'+str(I)]] = float(v_['AVP'])
            pb.x0[ix_['U'+str(I)]] = float(v_['UMAX'])
            pb.x0[ix_['A'+str(I)]] = float(v_['AA'])
            pb.x0[ix_['B'+str(I)]] = float(v_['BB'])
            pb.x0[ix_['C'+str(I)]] = float(v_['CC'])
        pb.x0[ix_['TF']] = float(240.0)
        pb.x0[ix_['W'+str(int(v_['NI']))]] = float(v_['WSS'])
        pb.x0[ix_['P'+str(int(v_['NI']))]] = float(v_['PFINAL'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en3S', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        [it,iet_,_] = s2mpj_ii( 'en3D2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        [it,iet_,_] = s2mpj_ii( 'eDSQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'en3PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            ename = 'WA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3S')
            ielftype = arrset(ielftype, ie, iet_["en3S"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'A'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'WB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3D2')
            ielftype = arrset(ielftype, ie, iet_["en3D2"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'PA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3D2')
            ielftype = arrset(ielftype, ie, iet_["en3D2"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'PB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3S')
            ielftype = arrset(ielftype, ie, iet_["en3S"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'DD'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eDSQ')
            ielftype = arrset(ielftype, ie, iet_["eDSQ"])
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'CA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
            ielftype = arrset(ielftype, ie, iet_["en3PR"])
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'A'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'CB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en3PR')
            ielftype = arrset(ielftype, ie, iet_["en3PR"])
            vname = 'C'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            ig = ig_['EW'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['WA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/NI']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['WB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-Z/NI']))
            ig = ig_['EP'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/NI']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-Z/NI']))
            ig = ig_['EA'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['EB'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DD'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['EC'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DD'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['CB'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-ZZ']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 200.0
#    Solution
# LO SOLTN               ????
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
        #%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        lxs = len(xscale);
        for j in np.arange(0,min(sA2,pb.n,len(xscale))):
            if not xscale[j] is None and xscale[j] != 0.0 and xscale[j] != 1.0:
                for i in find(pbm.A[:,j],lambda x:x!=0):
                      pbm.A[i,j] = pbm.A[i,j]/xscale[j]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-MY-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,0.02)
        return pbm

    @staticmethod
    def en3S(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        WSSMV4 = pbm.efpar[0]-EV_[3]
        f_   = EV_[0]*EV_[1]*EV_[2]*WSSMV4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]*WSSMV4
            g_[1] = EV_[0]*EV_[2]*WSSMV4
            g_[2] = EV_[0]*EV_[1]*WSSMV4
            g_[3] = -EV_[0]*EV_[1]*EV_[2]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = EV_[2]*WSSMV4
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]*WSSMV4
                H_[2,0] = H_[0,2]
                H_[0,3] = -EV_[1]*EV_[2]
                H_[3,0] = H_[0,3]
                H_[1,2] = EV_[0]*WSSMV4
                H_[2,1] = H_[1,2]
                H_[1,3] = -EV_[0]*EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = -EV_[0]*EV_[1]
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3D2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,5))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[3,3] = U_[3,3]+1
        U_[3,4] = U_[3,4]-2
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        IV_[3] = U_[3:4,:].dot(EV_)
        f_   = IV_[0]*IV_[1]*IV_[2]*IV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]*IV_[2]*IV_[3]
            g_[1] = IV_[0]*IV_[2]*IV_[3]
            g_[2] = IV_[0]*IV_[1]*IV_[3]
            g_[3] = IV_[0]*IV_[1]*IV_[2]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = IV_[2]*IV_[3]
                H_[1,0] = H_[0,1]
                H_[0,2] = IV_[1]*IV_[3]
                H_[2,0] = H_[0,2]
                H_[0,3] = IV_[1]*IV_[2]
                H_[3,0] = H_[0,3]
                H_[1,2] = IV_[0]*IV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = IV_[0]*IV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = IV_[0]*IV_[1]
                H_[3,2] = H_[2,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eDSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+2.000000e-01
        U_[0,1] = U_[0,1]+2.000000e-01
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3PR(pbm,nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

