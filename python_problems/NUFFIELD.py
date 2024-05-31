from s2xlib import *
class  NUFFIELD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NUFFIELD
#    *********
# 
#    A problem from economics.
#    Maximize a 2-D integral representing consumer surplus subject to 
#    linear and quadratic constraints representing incentive compatibility
# 
#    Let v( . , . ) : R^2 -> R, Omega = [a,a+1] x [a,a+1], and
#    the corners A, B, C, D be as follows:
# 
#            (a+1,a+1)
#        A *-----* B
#          |     |
#          |     |
#        D *-----* C
#        (a,a)  
# 
#    The problem is to maximize
# 
#       (a+1) line integral_{AB U BC} v(w)dw 
#        - a line integral_{CD U DA} v(w)dw
#        - 3 volume integral_{Omega} v(w)dw
# 
#    subject to v being symmetric (i.e., v(x,y) = v(y,x))
#               v(a,a) = 0
#               nabla_w v(w) >= 0
#               < e, nabla_w v(w) > <= 1
#         and   nabla_ww v(w) positive definite
# 
#    this last constraint is guaranteed by ensuring that
# 
#               d^2 v/dx^2 >= 0
#               d^2 v/dy^2 >= 0
#               ( d^2 v/dx^2 )( d^2 v/dy^2 ) >= ( d^2 v/dxdy )^2
# 
#    Symmetry is ensured by only considering v(x,y) for x <= y
# 
#    Here v(x,y) is the consumer surplus. that is if the consumer values good 
#    1 at x pounds and good 2 at y pounds then they will have a utility 
#    equivalent to v(x,y) pounds after being faced with the optimal monopoly 
#    pricing strategy. (Apparently, from this we can infer what the optimal 
#    pricing strategy was... ).
# 
#    More background is available from
# 
#    "Optimal Selling Strategies: When to haggle, when to hold firm",
#      Riley and Zeckhauser. The Quarterly Journal of Economics, 1983, and
# 
#    "Multidimensional Incentive Compatibility and Mechanism Design", 
#      McAfee and McMillan. The Journal of Economic Theory, 1988.
# 
#    Source: John Thanassoulis <john.thanassoulis@nuffield.oxford.ac.uk>
# 
#    Standard finite-differences are used to ap[proximate derivatives, and 
#    1- and 2-D trapezoidal rules to approximate integrals
# 
#    SIF input: Nick Gould, February 2001
# 
#    classification = "LQR2-AN-V-V"
# 
#    The parameter a
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'NUFFIELD'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'NUFFIELD'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['A'] = float(5.0);  #  SIF file default value
        else:
            v_['A'] = float(args[0])
        if nargin<2:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[1])
#           Alternative values for the SIF file parameters:
# IE N                   20            $-PARAMETER
# IE N                   30            $-PARAMETER
# IE N                   40            $-PARAMETER
# IE N                   100           $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = 1.0/v_['RN']
        v_['1/H'] = v_['RN']
        v_['-1/H'] = -1.0*v_['1/H']
        v_['H**2'] = v_['H']*v_['H']
        v_['1/H**2'] = v_['1/H']*v_['1/H']
        v_['-2/H**2'] = -2.0*v_['1/H**2']
        v_['1/H**4'] = v_['1/H**2']*v_['1/H**2']
        v_['A+1'] = 1.0+v_['A']
        v_['-A-1'] = -1.0*v_['A+1']
        v_['C2'] = 3.0*v_['H']
        v_['C3'] = 0.5*v_['C2']
        v_['C4'] = v_['C3']+v_['A']
        v_['C1'] = v_['C3']+v_['-A-1']
        v_['C5'] = -1.0+v_['C3']
        v_['C5'] = 0.5*v_['C5']
        v_['C6'] = 0.5*v_['C3']
        v_['C6'] = v_['C6']+v_['-A-1']
        v_['C6'] = 0.5*v_['C6']
        v_['C1'] = v_['C1']*v_['H']
        v_['C2'] = v_['C2']*v_['H']
        v_['C3'] = v_['C3']*v_['H']
        v_['C4'] = v_['C4']*v_['H']
        v_['C5'] = v_['C5']*v_['H']
        v_['C6'] = v_['C6']*v_['H']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(I)+1):
                [iv,ix_,_] = s2x_ii('V'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['V'+str(int(v_['N']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['C1'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2x_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['C2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['V'+str(I)+','+str(I)]
            pbm.A[ig,iv] = float(v_['C3'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['V'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(v_['C4'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['0']))]
        pbm.A[ig,iv] = float(v_['C5'])+pbm.A[ig,iv]
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['C6'])+pbm.A[ig,iv]
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['0']),int(I)+1):
                [ig,ig_,_] = s2x_ii('VX'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'VX'+str(I)+','+str(J))
                iv = ix_['V'+str(int(v_['I+1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('VV'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'VV'+str(I)+','+str(J))
                iv = ix_['V'+str(int(v_['I+1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        for J in range(int(v_['0']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2x_ii('VX'+str(int(v_['N']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'VX'+str(int(v_['N']))+','+str(J))
            iv = ix_['V'+str(int(v_['N']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('VX'+str(int(v_['N']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'VX'+str(int(v_['N']))+','+str(J))
            iv = ix_['V'+str(int(v_['N-1']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(J))
            iv = ix_['V'+str(int(v_['N']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(J),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(J))
            iv = ix_['V'+str(int(v_['N-1']))+','+str(J)]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VX'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'VX'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VX'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'VX'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['0']),int(v_['I-1'])+1):
                v_['J+1'] = 1+J
                [ig,ig_,_] = s2x_ii('VY'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'VY'+str(I)+','+str(J))
                iv = ix_['V'+str(I)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('VV'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'VV'+str(I)+','+str(J))
                iv = ix_['V'+str(I)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('VY'+str(I)+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'VY'+str(I)+','+str(I))
            iv = ix_['V'+str(int(v_['I+1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(I)]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('VV'+str(I)+','+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'VV'+str(I)+','+str(I))
            iv = ix_['V'+str(int(v_['I+1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(I)]
            pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VY'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'VY'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VY'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'VY'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['1/H'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('VV'+str(int(v_['N']))+','+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'VV'+str(int(v_['N']))+','+str(int(v_['N'])))
        iv = ix_['V'+str(int(v_['N']))+','+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(v_['-1/H'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['0']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2x_ii('VXX'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'VXX'+str(I)+','+str(J))
                iv = ix_['V'+str(int(v_['I+1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/H**2'])+pbm.A[ig,iv]
                iv = ix_['V'+str(int(v_['I-1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                [ig,ig_,_] = s2x_ii('VYY'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'VYY'+str(I)+','+str(J))
                iv = ix_['V'+str(I)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/H**2'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('VXX'+str(I)+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'VXX'+str(I)+','+str(I))
            iv = ix_['V'+str(int(v_['I+1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(I)]
            pbm.A[ig,iv] = float(v_['-2/H**2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('VYY'+str(I)+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'VYY'+str(I)+','+str(I))
            iv = ix_['V'+str(int(v_['I+1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(I)]
            pbm.A[ig,iv] = float(v_['-2/H**2'])+pbm.A[ig,iv]
            iv = ix_['V'+str(I)+','+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['1/H**2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(I)+1):
                [ig,ig_,_] = s2x_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
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
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(I)+1):
                pbm.gconst = arrset(pbm.gconst,ig_['VV'+str(I)+','+str(J)],float(1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['V'+str(int(v_['0']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['V'+str(int(v_['0']))+','+str(int(v_['0']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eCONVEX', iet_)
        elftv = loaset(elftv,it,0,'VIP1J')
        elftv = loaset(elftv,it,1,'VIJP1')
        elftv = loaset(elftv,it,2,'VIJ')
        elftv = loaset(elftv,it,3,'VIM1J')
        elftv = loaset(elftv,it,4,'VIJM1')
        elftv = loaset(elftv,it,5,'VIPJP')
        elftv = loaset(elftv,it,6,'VIPJM')
        elftv = loaset(elftv,it,7,'VIMJM')
        elftv = loaset(elftv,it,8,'VIMJP')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(I)+1):
                v_['J+1'] = 1+J
                v_['J-1'] = -1+J
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eCONVEX')
                ielftype = arrset(ielftype, ie, iet_["eCONVEX"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'V'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIP1J')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIM1J')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIJP1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIJM1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIJ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I+1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIPJP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIMJM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I+1']))+','+str(int(v_['J-1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIPJM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(int(v_['I-1']))+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='VIMJP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            for J in range(int(v_['1']),int(I)+1):
                ig = ig_['C'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/H**4']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCONVEX(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,9))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-2
        U_[0,3] = U_[0,3]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-2
        U_[1,4] = U_[1,4]+1
        U_[2,5] = U_[2,5]+2.500000e-01
        U_[2,7] = U_[2,7]+2.500000e-01
        U_[2,8] = U_[2,8]-2.500000e-01
        U_[2,6] = U_[2,6]-2.500000e-01
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        f_   = IV_[0]*IV_[1]-IV_[2]*IV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_[2] = -2.0*IV_[2]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_[2,2] = -2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

