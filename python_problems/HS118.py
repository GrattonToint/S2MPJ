from s2mpjlib import *
class  HS118(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS118
#    *********
# 
#    Source: problem 118 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: B Baudson, Jan 1990.
# 
#    classification = "QLR2-AN-15-17"
# 
#    Other useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS118'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['4'] = 4
        v_['8'] = 8
        v_['12'] = 12
        v_['15'] = 15
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['15'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for K in range(int(v_['0']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(int(v_['3K+1']))]
            self.A[ig,iv] = float(2.3)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K+2']))]
            self.A[ig,iv] = float(1.7)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K+3']))]
            self.A[ig,iv] = float(2.2)+self.A[ig,iv]
        for K in range(int(v_['1']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            v_['3K-2'] = -2+v_['3K']
            v_['3K-1'] = -1+v_['3K']
            [ig,ig_,_] = s2mpj_ii('A'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'A'+str(K))
            iv = ix_['X'+str(int(v_['3K+1']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K-2']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('B'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'B'+str(K))
            iv = ix_['X'+str(int(v_['3K+3']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(K),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(K))
            iv = ix_['X'+str(int(v_['3K+2']))]
            self.A[ig,iv] = float(1.0)+self.A[ig,iv]
            iv = ix_['X'+str(int(v_['3K-1']))]
            self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('D1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'D1')
        iv = ix_['X1']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X2']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X3']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('D2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'D2')
        iv = ix_['X4']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X5']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X6']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('D3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'D3')
        iv = ix_['X7']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X8']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X9']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('D4',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'D4')
        iv = ix_['X10']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X11']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X12']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('D5',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'D5')
        iv = ix_['X13']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X14']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        iv = ix_['X15']
        self.A[ig,iv] = float(1.0)+self.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames= cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for K in range(int(v_['1']),int(v_['4'])+1):
            self.gconst = arrset(self.gconst,ig_['A'+str(K)],float(-7.0))
            self.gconst = arrset(self.gconst,ig_['B'+str(K)],float(-7.0))
            self.gconst = arrset(self.gconst,ig_['C'+str(K)],float(-7.0))
        self.gconst = arrset(self.gconst,ig_['D1'],float(60.0))
        self.gconst = arrset(self.gconst,ig_['D2'],float(50.0))
        self.gconst = arrset(self.gconst,ig_['D3'],float(70.0))
        self.gconst = arrset(self.gconst,ig_['D4'],float(85.0))
        self.gconst = arrset(self.gconst,ig_['D5'],float(100.0))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        for K in range(int(v_['1']),int(v_['4'])+1):
            grange = arrset(grange,ig_['A'+str(K)],float(13.0))
            grange = arrset(grange,ig_['B'+str(K)],float(13.0))
            grange = arrset(grange,ig_['C'+str(K)],float(14.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.0)
        self.xupper = np.full((self.n,1),+float('inf'))
        self.xlower[ix_['X1']] = 8.0
        self.xupper[ix_['X1']] = 21.0
        self.xlower[ix_['X2']] = 43.0
        self.xupper[ix_['X2']] = 57.0
        self.xlower[ix_['X3']] = 3.0
        self.xupper[ix_['X3']] = 16.0
        for K in range(int(v_['1']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            self.xupper[ix_['X'+str(int(v_['3K+1']))]] = 90.0
            self.xupper[ix_['X'+str(int(v_['3K+2']))]] = 120.0
            self.xupper[ix_['X'+str(int(v_['3K+3']))]] = 60.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(20.0))
        if('X2' in ix_):
            self.x0[ix_['X2']] = float(55.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X2']),float(55.0)))
        if('X3' in ix_):
            self.x0[ix_['X3']] = float(15.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X3']),float(15.0)))
        if('X5' in ix_):
            self.x0[ix_['X5']] = float(60.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X5']),float(60.0)))
        if('X8' in ix_):
            self.x0[ix_['X8']] = float(60.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X8']),float(60.0)))
        if('X11' in ix_):
            self.x0[ix_['X11']] = float(60.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X11']),float(60.0)))
        if('X14' in ix_):
            self.x0[ix_['X14']] = float(60.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X14']),float(60.0)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['15'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,0.0,None,20.0)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['0']),int(v_['4'])+1):
            v_['3K'] = 3*K
            v_['3K+1'] = 1+v_['3K']
            v_['3K+2'] = 2+v_['3K']
            v_['3K+3'] = 3+v_['3K']
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3K+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(0.0001))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3K+2']))])
            self.grelw = loaset(self.grelw,ig,posel,float(0.0001))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3K+3']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(0.00015))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               664.82045
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nge)] = grange[gegrps]
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "QLR2-AN-15-17"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

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

