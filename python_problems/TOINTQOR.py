from s2mpjlib import *
class  TOINTQOR(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TOINTQOR
#    *********
# 
#    Toint's  Quadratic Operations Research problem
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
#    classification = "C-CQUR2-MN-50-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 9 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TOINTQOR'

    def __init__(self, *args): 
        import numpy as np
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('GA'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            v_['SCALE'] = 1.0/v_['ALPH'+str(I)]
            self.gscale = arrset(self.gscale,ig,float(v_['SCALE']))
        [ig,ig_,_] = s2mpj_ii('GB1',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X31']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB2',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X1']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB3',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X2']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB4',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X4']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X7']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB5',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X6']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB6',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X8']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB7',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X10']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB8',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X12']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB9',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X11']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X13']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X16']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X17']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB10',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X16']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X18']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X19']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB11',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X9']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X18']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X20']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB12',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X5']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X20']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X21']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB13',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X19']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X22']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X23']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X24']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB14',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X23']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X25']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X26']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB15',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X7']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X25']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X27']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X28']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB16',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X28']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X29']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X30']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB17',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X29']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X31']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X32']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB18',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X32']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X33']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X34']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB19',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X3']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X33']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X35']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB20',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X35']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X21']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X36']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB21',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X36']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X37']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X38']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB22',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X30']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X37']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X39']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB23',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X38']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X39']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X40']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB24',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X40']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X41']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X42']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB25',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X41']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X43']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X44']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X50']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB26',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X44']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X45']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X46']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X47']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB27',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X46']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X48']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB28',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X42']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X45']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X48']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X50']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X49']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB29',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X26']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X34']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X43']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB30',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X15']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X17']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X24']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        iv = ix_['X47']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB31',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X49']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB32',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X22']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GB33',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['X27']
        self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['33'])+1):
            v_['SCALE'] = 1.0/v_['BETA'+str(I)]
            [ig,ig_,_] = s2mpj_ii('GB'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            self.gscale = arrset(self.gscale,ig,float(v_['SCALE']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['33'])+1):
            self.gconst = arrset(self.gconst,ig_['GB'+str(I)],float(v_['D'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              1175.4722221
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "C-CQUR2-MN-50-0"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

