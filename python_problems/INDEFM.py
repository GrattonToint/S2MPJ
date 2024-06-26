from s2mpjlib import *
class  INDEFM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : INDEFM
#    *********
# 
#    Variant of INDEF, a nonconvex problem which has an indefinite Hessian 
#    at the starting point, by Luksan et al
# 
#    Source: problem 37 in
#    L. Luksan, C. Matonoha and J. Vlcek  
#    Modified CUTE problems for sparse unconstraoined optimization
#    Technical Report 1081
#    Institute of Computer Science
#    Academy of Science of the Czech Republic
# 
#    based on the original problem by N. Gould
# 
#    SIF input: Nick Gould, June, 2013
# 
#    classification = "OUR2-AN-V-0"
# 
#    The number of variables is N.
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'INDEFM'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# IE N                   100000         $-PARAMETER
        if nargin<2:
            v_['ALPHA'] = float(0.5);  #  SIF file default value
        else:
            v_['ALPHA'] = float(args[1])
# RE ALPHA               1.0            $-PARAMETER
# RE ALPHA               10.0           $-PARAMETER
# RE ALPHA               100.0          $-PARAMETER
# RE ALPHA               1000.0         $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['N-1'] = -1+v_['N']
        v_['N+1'] = 1+v_['N']
        v_['RN+1'] = float(v_['N+1'])
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
            [ig,ig_,_] = s2mpj_ii('SIN'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            [ig,ig_,_] = s2mpj_ii('COS'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            self.A[ig,iv] = float(2.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['N']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['1']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['T'] = v_['RI']/v_['RN+1']
            self.x0[ix_['X'+str(I)]] = float(v_['T'])
        pass
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gCOS',igt_)
        [it,igt_,_] = s2mpj_ii('gCOS',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'ALPHA')
        [it,igt_,_] = s2mpj_ii('gSIN',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        self.grpar   = []
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            ig = ig_['COS'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gCOS')
            posgp = np.where(grftp[igt_[self.grftype[ig]]]=='ALPHA')[0]
            self.grpar =loaset(self.grpar,ig,posgp[0],float(v_['ALPHA']))
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['SIN'+str(I)]
            self.grftype = arrset(self.grftype,ig,'gSIN')
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               ??
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "OUR2-AN-V-0"
# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gCOS(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= self.grpar[igr_][0]*np.cos(GVAR_)
        if nargout>1:
            g_ = -self.grpar[igr_][0]*np.sin(GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -self.grpar[igr_][0]*np.cos(GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def gSIN(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 100.0*np.sin(0.01*GVAR_)
        if nargout>1:
            g_ = np.cos(0.01*GVAR_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.01*np.sin(0.01*GVAR_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

