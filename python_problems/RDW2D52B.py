from s2mpjlib import *
class  RDW2D52B(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RDW2D52B
#    *********
# 
#    A finite-element approximation to the distributed optimal control problem
# 
#       min 1/2||u-v||_L2^2 + beta ||f||_L2^2
# 
#    subject to - nabla^2 u = f
# 
#    where v is given on and within the boundary of a unit [0,1] box in 
#    2 dimensions, and u = v on its boundary. The discretization uses 
#    quadrilateral elememts. There are simple bounds on both the controls 
#    f and states u
# 
#    The problem is stated as a quadratic program
# 
#    Source:  example 5.2 in 
#     T. Rees, H. S. Dollar and A. J. Wathen
#     "Optimal solvers for PDE-constrained optimization"
#     SIAM J. Sci. Comp. (to appear) 2009
# 
#    with the control bounds as specified in 
# 
#     M. Stoll and A. J. Wathen
#     "Preconditioning for PDE constrained optimization with 
#      control constraints"
#     OUCL Technical Report 2009
# 
#    SIF input: Nick Gould, May 2009
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "QLR2-AN-V-V"
# 
#    Number of nodes in each direction (a power of 2)
# 
#           Alternative values for the SIF file parameters:
# IE N                   2             $-PARAMETER
# IE N                   4             $-PARAMETER
# IE N                   8             $-PARAMETER
# IE N                   16            $-PARAMETER
# IE N                   32            $-PARAMETER
# IE N                   64            $-PARAMETER
# IE N                   128           $-PARAMETER
# IE N                   256           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RDW2D52B'

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
            v_['N'] = int(4);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   512           $-PARAMETER
# IE N                   1024          $-PARAMETER
# IE N                   2048          $-PARAMETER
# IE N                   4096          $-PARAMETER
# IE N                   8192          $-PARAMETER
# IE N                   16384         $-PARAMETER
        if nargin<2:
            v_['BETA'] = float(0.005);  #  SIF file default value
        else:
            v_['BETA'] = float(args[1])
        v_['ZERO'] = 0.0
        v_['ONE'] = 1.0
        v_['TWO'] = 2.0
        v_['SIX'] = 6.0
        v_['THIRTYSIX'] = 36.0
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N-'] = -1+v_['N']
        v_['N-1'] = -1+v_['N']
        v_['N-2'] = -2+v_['N']
        v_['N/2'] = int(np.fix(v_['N']/v_['2']))
        v_['N/2+1'] = v_['N/2']+v_['1']
        v_['RN'] = float(v_['N'])
        v_['H'] = v_['ONE']/v_['RN']
        v_['H**2'] = v_['H']*v_['H']
        v_['H**2/36'] = v_['H**2']/v_['THIRTYSIX']
        v_['-H**2/36'] = -1.0*v_['H**2/36']
        v_['2BETA'] = 2.0*v_['BETA']
        v_['2BH**2/36'] = v_['2BETA']*v_['H**2/36']
        v_['1/6'] = v_['ONE']/v_['SIX']
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['X'] = v_['RI']*v_['H']
            v_['X-'] = -0.5+v_['X']
            v_['X-**2'] = v_['X-']*v_['X-']
            for J in range(int(v_['0']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['Y'] = v_['RJ']*v_['H']
                v_['Y-'] = -0.5+v_['Y']
                v_['Y-**2'] = v_['Y-']*v_['Y-']
                v_['SS'] = v_['X-**2']+v_['Y-**2']
                v_['-64SS'] = -64.0*v_['SS']
                v_['V'] = np.exp(v_['-64SS'])
                v_['V'+str(I)+','+str(J)] = v_['V']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('F'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'F'+str(I)+','+str(J))
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['N-1'])+1):
                [ig,ig_,_] = s2mpj_ii('L'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'L'+str(I)+','+str(J))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['U'+str(int(v_['0']))+','+str(int(v_['0']))]] = (v_['V'+
             str(int(v_['0']))+','+str(int(v_['0']))])
        pb.xupper[ix_['U'+str(int(v_['0']))+','+str(int(v_['0']))]] = (v_['V'+
             str(int(v_['0']))+','+str(int(v_['0']))])
        pb.xlower[ix_['U'+str(int(v_['N']))+','+str(int(v_['0']))]] = (v_['V'+
             str(int(v_['N']))+','+str(int(v_['0']))])
        pb.xupper[ix_['U'+str(int(v_['N']))+','+str(int(v_['0']))]] = (v_['V'+
             str(int(v_['N']))+','+str(int(v_['0']))])
        pb.xlower[ix_['U'+str(int(v_['0']))+','+str(int(v_['N']))]] = (v_['V'+
             str(int(v_['0']))+','+str(int(v_['N']))])
        pb.xupper[ix_['U'+str(int(v_['0']))+','+str(int(v_['N']))]] = (v_['V'+
             str(int(v_['0']))+','+str(int(v_['N']))])
        pb.xlower[ix_['U'+str(int(v_['N']))+','+str(int(v_['N']))]] = (v_['V'+
             str(int(v_['N']))+','+str(int(v_['N']))])
        pb.xupper[ix_['U'+str(int(v_['N']))+','+str(int(v_['N']))]] = (v_['V'+
             str(int(v_['N']))+','+str(int(v_['N']))])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pb.xlower[ix_['U'+str(int(v_['0']))+','+str(I)]] = (v_['V'+str(int(v_['0']))+
                 ','+str(I)])
            pb.xupper[ix_['U'+str(int(v_['0']))+','+str(I)]] = (v_['V'+str(int(v_['0']))+
                 ','+str(I)])
            pb.xlower[ix_['U'+str(int(v_['N']))+','+str(I)]] = (v_['V'+str(int(v_['N']))+
                 ','+str(I)])
            pb.xupper[ix_['U'+str(int(v_['N']))+','+str(I)]] = (v_['V'+str(int(v_['N']))+
                 ','+str(I)])
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['0']))]] = (v_['V'+str(I)+
                 ','+str(int(v_['0']))])
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['0']))]] = (v_['V'+str(I)+
                 ','+str(int(v_['0']))])
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['N']))]] = (v_['V'+str(I)+
                 ','+str(int(v_['N']))])
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['N']))]] = (v_['V'+str(I)+
                 ','+str(int(v_['N']))])
        pb.xlower[ix_['F'+str(int(v_['0']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['F'+str(int(v_['0']))+','+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['F'+str(int(v_['N']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['F'+str(int(v_['N']))+','+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['F'+str(int(v_['0']))+','+str(int(v_['N']))]] = 0.0
        pb.xupper[ix_['F'+str(int(v_['0']))+','+str(int(v_['N']))]] = 0.0
        pb.xlower[ix_['F'+str(int(v_['N']))+','+str(int(v_['N']))]] = 0.0
        pb.xupper[ix_['F'+str(int(v_['N']))+','+str(int(v_['N']))]] = 0.0
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            pb.xlower[ix_['F'+str(int(v_['0']))+','+str(I)]] = 0.0
            pb.xupper[ix_['F'+str(int(v_['0']))+','+str(I)]] = 0.0
            pb.xlower[ix_['F'+str(int(v_['N']))+','+str(I)]] = 0.0
            pb.xupper[ix_['F'+str(int(v_['N']))+','+str(I)]] = 0.0
            pb.xlower[ix_['F'+str(I)+','+str(int(v_['0']))]] = 0.0
            pb.xupper[ix_['F'+str(I)+','+str(int(v_['0']))]] = 0.0
            pb.xlower[ix_['F'+str(I)+','+str(int(v_['N']))]] = 0.0
            pb.xupper[ix_['F'+str(I)+','+str(int(v_['N']))]] = 0.0
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(v_['N-1'])+1):
                pb.xupper[ix_['U'+str(I)+','+str(J)]] = 0.01
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RI'] = float(I)
            v_['X1'] = v_['RI']*v_['H']
            v_['X1**2'] = v_['X1']*v_['X1']
            v_['-X1**2'] = -1.0*v_['X1**2']
            v_['2-X1'] = v_['TWO']-v_['X1']
            v_['0.1(2-X1)'] = 0.1*v_['2-X1']
            for J in range(int(v_['1']),int(v_['N-1'])+1):
                v_['RJ'] = float(J)
                v_['X2'] = v_['RJ']*v_['H']
                v_['X2**2'] = v_['X2']*v_['X2']
                v_['ARG'] = v_['-X1**2']-v_['X2**2']
                v_['EARG'] = np.exp(v_['ARG'])
                v_['UA'] = v_['0.1(2-X1)']*v_['EARG']
                pb.xlower[ix_['F'+str(I)+','+str(J)]] = v_['UA']
            for J in range(int(v_['1']),int(v_['N/2'])+1):
                pb.xupper[ix_['F'+str(I)+','+str(J)]] = 0.6
            for J in range(int(v_['N/2+1']),int(v_['N-1'])+1):
                pb.xupper[ix_['F'+str(I)+','+str(J)]] = 0.9
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                v_['I+J'] = I+J
                v_['RI+J'] = float(v_['I+J'])
                v_['RI+J/N'] = v_['RI+J']/v_['RN']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eM', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'U4')
        elftp = []
        elftp = loaset(elftp,it,0,'V1')
        elftp = loaset(elftp,it,1,'V2')
        elftp = loaset(elftp,it,2,'V3')
        elftp = loaset(elftp,it,3,'V4')
        [it,iet_,_] = s2mpj_ii( 'eM0', iet_)
        elftv = loaset(elftv,it,0,'F1')
        elftv = loaset(elftv,it,1,'F2')
        elftv = loaset(elftv,it,2,'F3')
        elftv = loaset(elftv,it,3,'F4')
        [it,iet_,_] = s2mpj_ii( 'eA', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'U4')
        [it,iet_,_] = s2mpj_ii( 'eB', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'U4')
        [it,iet_,_] = s2mpj_ii( 'eC', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'U4')
        [it,iet_,_] = s2mpj_ii( 'eD', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        elftv = loaset(elftv,it,3,'U4')
        [it,iet_,_] = s2mpj_ii( 'eP', iet_)
        elftv = loaset(elftv,it,0,'F1')
        elftv = loaset(elftv,it,1,'F2')
        elftv = loaset(elftv,it,2,'F3')
        elftv = loaset(elftv,it,3,'F4')
        [it,iet_,_] = s2mpj_ii( 'eQ', iet_)
        elftv = loaset(elftv,it,0,'F1')
        elftv = loaset(elftv,it,1,'F2')
        elftv = loaset(elftv,it,2,'F3')
        elftv = loaset(elftv,it,3,'F4')
        [it,iet_,_] = s2mpj_ii( 'eR', iet_)
        elftv = loaset(elftv,it,0,'F1')
        elftv = loaset(elftv,it,1,'F2')
        elftv = loaset(elftv,it,2,'F3')
        elftv = loaset(elftv,it,3,'F4')
        [it,iet_,_] = s2mpj_ii( 'eS', iet_)
        elftv = loaset(elftv,it,0,'F1')
        elftv = loaset(elftv,it,1,'F2')
        elftv = loaset(elftv,it,2,'F3')
        elftv = loaset(elftv,it,3,'F4')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+'] = I+v_['1']
            for J in range(int(v_['0']),int(v_['N-1'])+1):
                v_['J+'] = J+v_['1']
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eM')
                ielftype = arrset(ielftype, ie, iet_["eM"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='V1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['V'+str(I)+','+str(J)]))
                posep = find(elftp[ielftype[ie]],lambda x:x=='V2')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['V'+str(I)+','+str(int(v_['J+']))])))
                posep = find(elftp[ielftype[ie]],lambda x:x=='V3')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['V'+str(int(v_['I+']))+','+str(J)])))
                posep = find(elftp[ielftype[ie]],lambda x:x=='V4')
                pbm.elpar  = (
                      loaset(pbm.elpar,ie,posep[0],float(v_['V'+str(int(v_['I+']))+','+str(int(v_['J+']))])))
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+'] = I+v_['1']
            for J in range(int(v_['0']),int(v_['N-1'])+1):
                v_['J+'] = J+v_['1']
                ename = 'F'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eM0')
                ielftype = arrset(ielftype, ie, iet_["eM0"])
                vname = 'F'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+'] = I+v_['1']
            for J in range(int(v_['0']),int(v_['N-1'])+1):
                v_['J+'] = J+v_['1']
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eA')
                ielftype = arrset(ielftype, ie, iet_["eA"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'B'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eB')
                ielftype = arrset(ielftype, ie, iet_["eB"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eC')
                ielftype = arrset(ielftype, ie, iet_["eC"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'D'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eD')
                ielftype = arrset(ielftype, ie, iet_["eD"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+'] = I+v_['1']
            for J in range(int(v_['0']),int(v_['N-1'])+1):
                v_['J+'] = J+v_['1']
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eP')
                ielftype = arrset(ielftype, ie, iet_["eP"])
                vname = 'F'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Q'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eQ')
                ielftype = arrset(ielftype, ie, iet_["eQ"])
                vname = 'F'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'R'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eR')
                ielftype = arrset(ielftype, ie, iet_["eR"])
                vname = 'F'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eS')
                ielftype = arrset(ielftype, ie, iet_["eS"])
                vname = 'F'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'F'+str(int(v_['I+']))+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='F4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            for J in range(int(v_['0']),int(v_['N-1'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['H**2/36']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['2BH**2/36']))
        for I in range(int(v_['1']),int(v_['N-2'])+1):
            v_['I+'] = I+v_['1']
            for J in range(int(v_['1']),int(v_['N-2'])+1):
                v_['J+'] = J+v_['1']
                ig = ig_['L'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
                ig = ig_['L'+str(I)+','+str(int(v_['J+']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
                ig = ig_['L'+str(int(v_['I+']))+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
                ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['J+']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
                ig = ig_['L'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
                ig = ig_['L'+str(I)+','+str(int(v_['J+']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
                ig = ig_['L'+str(int(v_['I+']))+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['R'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
                ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['J+']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        for I in range(int(v_['1']),int(v_['N-2'])+1):
            v_['I+'] = I+v_['1']
            ig = ig_['L'+str(I)+','+str(int(v_['N-']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(int(v_['N-']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(I)+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(int(v_['0']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['N-']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(int(v_['N-']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)+','+str(int(v_['0']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(I)+','+str(int(v_['N-']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(int(v_['N-']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(I)+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(int(v_['0']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['N-']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['R'+str(I)+','+str(int(v_['N-']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(int(v_['I+']))+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(int(v_['0']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        for J in range(int(v_['1']),int(v_['N-2'])+1):
            v_['J+'] = J+v_['1']
            ig = ig_['L'+str(int(v_['N-']))+','+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['N-']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['J+']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['B'+str(int(v_['N-']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['1']))+','+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(int(v_['0']))+','+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['J+']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['0']))+','+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
            ig = ig_['L'+str(int(v_['N-']))+','+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['P'+str(int(v_['N-']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['J+']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['Q'+str(int(v_['N-']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(int(v_['1']))+','+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['R'+str(int(v_['0']))+','+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
            ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['J+']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['0']))+','+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['N-']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['A'+str(int(v_['N-']))+','+str(int(v_['N-']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
        ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['B'+str(int(v_['N-']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
        ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['N-']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['C'+str(int(v_['0']))+','+str(int(v_['N-']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
        ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['0']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/6']))
        ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['N-']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['P'+str(int(v_['N-']))+','+str(int(v_['N-']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        ig = ig_['L'+str(int(v_['N-']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['Q'+str(int(v_['N-']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['N-']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['R'+str(int(v_['0']))+','+str(int(v_['N-']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        ig = ig_['L'+str(int(v_['1']))+','+str(int(v_['1']))]
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['0']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-H**2/36']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eM(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        UV1 = EV_[0]-pbm.elpar[iel_][0]
        UV2 = EV_[1]-pbm.elpar[iel_][1]
        UV3 = EV_[2]-pbm.elpar[iel_][2]
        UV4 = EV_[3]-pbm.elpar[iel_][3]
        f_   = (2.0*UV1**2+2.0*UV2**2+2.0*UV3**2+2.0*UV4**2+2.0*UV1*UV2+2.0*UV1*UV3+
             UV1*UV4+UV2*UV3+2.0*UV2*UV4+2.0*UV3*UV4)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*UV1+2.0*UV2+2.0*UV3+UV4
            g_[1] = 2.0*UV1+4.0*UV2+UV3+2.0*UV4
            g_[2] = 2.0*UV1+UV2+4.0*UV3+2.0*UV4
            g_[3] = UV1+2.0*UV2+2.0*UV3+4.0*UV4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 4.0
                H_[0,1] = 2.0
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0
                H_[2,0] = H_[0,2]
                H_[0,3] = 1.0
                H_[3,0] = H_[0,3]
                H_[1,1] = 4.0
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0
                H_[3,1] = H_[1,3]
                H_[2,2] = 4.0
                H_[2,3] = 2.0
                H_[3,2] = H_[2,3]
                H_[3,3] = 4.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eM0(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_    = (
              2.0*EV_[0]**2+2.0*EV_[1]**2+2.0*EV_[2]**2+2.0*EV_[3]**2+2.0*EV_[0]*EV_[1]+2.0*EV_[0]*EV_[2]+EV_[0]*EV_[3]+EV_[1]*EV_[2]+2.0*EV_[1]*EV_[3]+2.0*EV_[2]*EV_[3])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*EV_[0]+2.0*EV_[1]+2.0*EV_[2]+EV_[3]
            g_[1] = 2.0*EV_[0]+4.0*EV_[1]+EV_[2]+2.0*EV_[3]
            g_[2] = 2.0*EV_[0]+EV_[1]+4.0*EV_[2]+2.0*EV_[3]
            g_[3] = EV_[0]+2.0*EV_[1]+2.0*EV_[2]+4.0*EV_[3]
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 4.0
                H_[0,1] = 2.0
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0
                H_[2,0] = H_[0,2]
                H_[0,3] = 1.0
                H_[3,0] = H_[0,3]
                H_[1,1] = 4.0
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0
                H_[3,1] = H_[1,3]
                H_[2,2] = 4.0
                H_[2,3] = 2.0
                H_[3,2] = H_[2,3]
                H_[3,3] = 4.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = 4.0
        C2 = -1.0
        C3 = -1.0
        C4 = -2.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = -1.0
        C2 = 4.0
        C3 = -2.0
        C4 = -1.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eC(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = -1.0
        C2 = -2.0
        C3 = 4.0
        C4 = -1.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = -2.0
        C2 = -1.0
        C3 = -1.0
        C4 = 4.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = 4.0
        C2 = 2.0
        C3 = 2.0
        C4 = 1.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = 2.0
        C2 = 4.0
        C3 = 1.0
        C4 = 2.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = 2.0
        C2 = 1.0
        C3 = 4.0
        C4 = 2.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eS(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        C1 = 1.0
        C2 = 2.0
        C3 = 2.0
        C4 = 4.0
        f_   = C1*EV_[0]+C2*EV_[1]+C3*EV_[2]+C4*EV_[3]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1
            g_[1] = C2
            g_[2] = C3
            g_[3] = C4
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = 0.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

