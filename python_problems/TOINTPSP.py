from s2mpjlib import *
class  TOINTPSP(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TOINTPSP
#    *********
# 
#    Toint's PSP Operations Research problem
# 
#    Source:
#    Ph.L. Toint,
#    "Some numerical results using a sparse matrix updating formula in
#    unconstrained optimization",
#    Mathematics of Computation 32(1):839-852, 1978.
# 
#    See also Buckley#55 (p.94) (With a slightly lower optimal value?)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "OUR2-AN-50-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TOINTPSP'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'TOINTPSP'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 50
        v_['ALPH1'] = 1.25
        v_['ALPH2'] = 1.40
        v_['ALPH3'] = 2.40
        v_['ALPH4'] = 1.40
        v_['ALPH5'] = 1.75
        v_['ALPH6'] = 1.20
        v_['ALPH7'] = 2.25
        v_['ALPH8'] = 1.20
        v_['ALPH9'] = 1.00
        v_['ALPH10'] = 1.10
        v_['ALPH11'] = 1.50
        v_['ALPH12'] = 1.60
        v_['ALPH13'] = 1.25
        v_['ALPH14'] = 1.25
        v_['ALPH15'] = 1.20
        v_['ALPH16'] = 1.20
        v_['ALPH17'] = 1.40
        v_['ALPH18'] = 0.50
        v_['ALPH19'] = 0.50
        v_['ALPH20'] = 1.25
        v_['ALPH21'] = 1.80
        v_['ALPH22'] = 0.75
        v_['ALPH23'] = 1.25
        v_['ALPH24'] = 1.40
        v_['ALPH25'] = 1.60
        v_['ALPH26'] = 2.00
        v_['ALPH27'] = 1.00
        v_['ALPH28'] = 1.60
        v_['ALPH29'] = 1.25
        v_['ALPH30'] = 2.75
        v_['ALPH31'] = 1.25
        v_['ALPH32'] = 1.25
        v_['ALPH33'] = 1.25
        v_['ALPH34'] = 3.00
        v_['ALPH35'] = 1.50
        v_['ALPH36'] = 2.00
        v_['ALPH37'] = 1.25
        v_['ALPH38'] = 1.40
        v_['ALPH39'] = 1.80
        v_['ALPH40'] = 1.50
        v_['ALPH41'] = 2.20
        v_['ALPH42'] = 1.40
        v_['ALPH43'] = 1.50
        v_['ALPH44'] = 1.25
        v_['ALPH45'] = 2.00
        v_['ALPH46'] = 1.50
        v_['ALPH47'] = 1.25
        v_['ALPH48'] = 1.40
        v_['ALPH49'] = 0.60
        v_['ALPH50'] = 1.50
        v_['BETA1'] = 1.0
        v_['BETA2'] = 1.5
        v_['BETA3'] = 1.0
        v_['BETA4'] = 0.1
        v_['BETA5'] = 1.5
        v_['BETA6'] = 2.0
        v_['BETA7'] = 1.0
        v_['BETA8'] = 1.5
        v_['BETA9'] = 3.0
        v_['BETA10'] = 2.0
        v_['BETA11'] = 1.0
        v_['BETA12'] = 3.0
        v_['BETA13'] = 0.1
        v_['BETA14'] = 1.5
        v_['BETA15'] = 0.15
        v_['BETA16'] = 2.0
        v_['BETA17'] = 1.0
        v_['BETA18'] = 0.1
        v_['BETA19'] = 3.0
        v_['BETA20'] = 0.1
        v_['BETA21'] = 1.2
        v_['BETA22'] = 1.0
        v_['BETA23'] = 0.1
        v_['BETA24'] = 2.0
        v_['BETA25'] = 1.2
        v_['BETA26'] = 3.0
        v_['BETA27'] = 1.5
        v_['BETA28'] = 3.0
        v_['BETA29'] = 2.0
        v_['BETA30'] = 1.0
        v_['BETA31'] = 1.2
        v_['BETA32'] = 2.0
        v_['BETA33'] = 1.0
        v_['D1'] = -5.0
        v_['D2'] = -5.0
        v_['D3'] = -5.0
        v_['D4'] = -2.5
        v_['D5'] = -6.0
        v_['D6'] = -6.0
        v_['D7'] = -5.0
        v_['D8'] = -6.0
        v_['D9'] = -10.0
        v_['D10'] = -6.0
        v_['D11'] = -5.0
        v_['D12'] = -9.0
        v_['D13'] = -2.0
        v_['D14'] = -7.0
        v_['D15'] = -2.5
        v_['D16'] = -6.0
        v_['D17'] = -5.0
        v_['D18'] = -2.0
        v_['D19'] = -9.0
        v_['D20'] = -2.0
        v_['D21'] = -5.0
        v_['D22'] = -5.0
        v_['D23'] = -2.5
        v_['D24'] = -5.0
        v_['D25'] = -6.0
        v_['D26'] = -10.0
        v_['D27'] = -7.0
        v_['D28'] = -10.0
        v_['D29'] = -6.0
        v_['D30'] = -5.0
        v_['D31'] = -4.0
        v_['D32'] = -4.0
        v_['D33'] = -4.0
        v_['1'] = 1
        v_['33'] = 33
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('GA'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            v_['SCALE'] = 1.0/v_['ALPH'+str(I)]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['SCALE']))
        [ig,ig_,_] = s2mpj_ii('GB1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X31']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X1']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X2']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X5']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X4']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X6']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X7']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X6']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X8']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X9']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB6',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X8']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X10']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X11']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB7',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X10']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X12']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB8',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X12']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X15']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB9',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X11']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X13']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X14']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X16']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB10',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X16']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X19']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB11',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X9']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X18']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB12',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X5']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X20']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X21']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB13',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X19']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X22']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X23']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X24']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB14',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X23']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X25']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X26']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB15',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X7']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X25']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X27']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X28']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB16',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X28']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X29']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X30']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB17',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X29']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X31']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X32']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB18',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X32']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X33']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X34']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB19',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X33']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X35']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB20',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X35']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X21']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X36']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB21',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X36']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X37']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X38']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB22',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X30']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X37']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X39']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB23',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X38']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X39']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X40']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB24',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X40']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X41']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X42']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB25',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X41']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X43']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X44']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X50']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB26',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X44']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X45']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X46']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['X47']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB27',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X46']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X48']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB28',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X42']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X45']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X48']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X50']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X49']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB29',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X26']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X34']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X43']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB30',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X15']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X17']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X24']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        iv = ix_['X47']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB31',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X49']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB32',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X22']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB33',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X27']
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['33'])+1):
            v_['SCALE'] = 1.0/v_['BETA'+str(I)]
            [ig,ig_,_] = s2mpj_ii('GB'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['SCALE']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['33'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['GB'+str(I)],float(v_['D'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gACT',igt_)
        [it,igt_,_] = s2mpj_ii('gBBT',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['GA'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gACT')
        for I in range(int(v_['1']),int(v_['33'])+1):
            ig = ig_['GB'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gBBT')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              225.56040942
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-50-0"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gACT(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= (GVAR_-5.0)**2
        if nargout>1:
            g_ = 2.0*GVAR_-10.0
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gBBT(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        TPOS = GVAR_>=0.1
        if TPOS!=0:
            FF = 1.0/GVAR_
        if TPOS==0:
            FF = 20.0-100.0*GVAR_
        if TPOS!=0:
            GG = -1.0/GVAR_**2
        if TPOS==0:
            GG = -100.0
        if TPOS!=0:
            HH = 2.0/GVAR_**3
        if TPOS==0:
            HH = 0.0
        f_= FF
        if nargout>1:
            g_ = GG
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = HH
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

