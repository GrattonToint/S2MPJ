from s2mpjlib import *
class  LISWET7(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LISWET7
#    *********
# 
#    A k-convex approximation problem posed as a 
#    convex quadratic problem, with variable dimensions.
# 
#    Formulation:
#    -----------
# 
#                 n+k             2
#    minimize 1/2 sum ( x  - c  )
#                 i=1    i    i
# 
#    subject to
# 
#                  k              k-i
#                 sum ( k ) ( -1 )    x     > 0
#                 i=0 ( i )            j+i  = 
# 
#    where c  = g( t ) + small perturbation, t  = (i-1)/(n+k-1)
#           i       i                         i 
# 
#    Case 7: g(t) = sin(pi t)
# 
#    NB. Perturbations are not random as Li and Swetits's 
#        random number generator is undefined.
# 
#    Source:
#    W. Li and J. Swetits,
#    "A Newton method for convex regression, data smoothing and
#    quadratic programming with bounded constraints",
#    SIAM J. Optimization 3 (3) pp 466-488, 1993.
# 
#    SIF input: Nick Gould, August 1994.
# 
#    classification = "QLR2-AN-V-V"
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER 103 variables original value 
# IE K                   3              $-PARAMETER original value
# 
# IE N                   100            $-PARAMETER 104 variables    
# IE K                   4              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 105 variables    
# IE K                   5              $-PARAMETER
# 
# IE N                   100            $-PARAMETER 106 variables    
# IE K                   6              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 402 variables    
# IE K                   2              $-PARAMETER
# 
# IE N                   400            $-PARAMETER 403 variables    
# IE K                   3              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2001 variables    
# IE K                   1              $-PARAMETER
# 
# IE N                   2000           $-PARAMETER 2002 variables    
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LISWET7'

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
            v_['N'] = int(50);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE K                   2              $-PARAMETER
        if nargin<2:
            v_['K'] = int(3);  #  SIF file default value
        else:
            v_['K'] = int(args[1])
# IE N                   10000          $-PARAMETER 10001 variables    
# IE K                   1              $-PARAMETER
# IE N                   10000          $-PARAMETER 10002 variables    
# IE K                   2              $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['HALF'] = 0.5
        v_['N+K'] = v_['N']+v_['K']
        v_['N+K-1'] = -1+v_['N+K']
        v_['RN+K-1'] = float(v_['N+K-1'])
        v_['CONST'] = 0.0
        v_['B'+str(int(v_['0']))] = v_['ONE']
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['B'+str(I)] = v_['B'+str(int(v_['I-1']))]*v_['RI']
        v_['C'+str(int(v_['0']))] = v_['ONE']
        v_['PLUSMINUS'] = v_['ONE']
        for I in range(int(v_['1']),int(v_['K'])+1):
            v_['K-I'] = v_['K']-I
            v_['PLUSMINUS'] = -1.0*v_['PLUSMINUS']
            v_['C'+str(I)] = v_['B'+str(int(v_['K']))]/v_['B'+str(I)]
            v_['C'+str(I)] = v_['C'+str(I)]/v_['B'+str(int(v_['K-I']))]
            v_['C'+str(I)] = v_['C'+str(I)]*v_['PLUSMINUS']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI-1'] = float(v_['I-1'])
            v_['TI'] = v_['RI-1']/v_['RN+K-1']
            v_['PI/4'] = np.arctan(1.0)
            v_['PI'] = 4.0*v_['PI/4']
            v_['PIT'] = v_['PI']*v_['TI']
            v_['GT'] = np.sin(v_['PIT'])
            v_['RANDOM'] = np.sin(v_['RI'])
            v_['RANDOM'] = 0.1*v_['RANDOM']
            v_['CI'] = v_['GT']+v_['RANDOM']
            v_['-CI'] = -1.0*v_['CI']
            v_['-CI*CI'] = v_['-CI']*v_['CI']
            v_['CONST'] = v_['CONST']+v_['-CI*CI']
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['-CI'])+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['J+K'] = J+v_['K']
            for I in range(int(v_['0']),int(v_['K'])+1):
                v_['J+K-I'] = v_['J+K']-I
                [ig,ig_,_] = s2mpj_ii('CON'+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'CON'+str(J))
                iv = ix_['X'+str(int(v_['J+K-I']))]
                pbm.A[ig,iv] = float(v_['C'+str(I)])+pbm.A[ig,iv]
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
        v_['CONST'] = v_['HALF']*v_['CONST']
        pbm.gconst = arrset(pbm.gconst,ig_['OBJ'],float(v_['CONST']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            ename = 'XSQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            pb.x0 = np.zeros((pb.n,1))
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N+K'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['XSQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 5.0e-1*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

