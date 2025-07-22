from s2mpjlib import *
class  HS67(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS67
#    *********
# 
#    Source: problem 67 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 8 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn & Nick Gould, April 1991.
# 
#    classification = "C-COOI2-AN-3-14"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 22 VII 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS67'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['A'+str(int(v_['1']))] = 0.0
        v_['A'+str(int(v_['2']))] = 0.0
        v_['A'+str(int(v_['3']))] = 85.0
        v_['A'+str(int(v_['4']))] = 90.0
        v_['A'+str(int(v_['5']))] = 3.0
        v_['A'+str(int(v_['6']))] = 0.01
        v_['A'+str(int(v_['7']))] = 145.0
        v_['A'+str(int(v_['8']))] = 5000.0
        v_['A'+str(int(v_['9']))] = 2000.0
        v_['A'+str(int(v_['10']))] = 93.0
        v_['A'+str(int(v_['11']))] = 95.0
        v_['A'+str(int(v_['12']))] = 12.0
        v_['A'+str(int(v_['13']))] = 4.0
        v_['A'+str(int(v_['14']))] = 162.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['3'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(5.04))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['2']))]])
        valA = np.append(valA,float(0.035))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['3']))]])
        valA = np.append(valA,float(10.0))
        for I in range(int(v_['1']),int(v_['7'])+1):
            [ig,ig_,_] = s2mpj_ii('AG'+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'AG'+str(I))
            [ig,ig_,_] = s2mpj_ii('AL'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'AL'+str(I))
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
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['7'])+1):
            v_['I+7'] = 7+I
            self.gconst = arrset(self.gconst,ig_['AG'+str(I)],float(v_['A'+str(I)]))
            self.gconst  = (
                  arrset(self.gconst,ig_['AL'+str(I)],float(v_['A'+str(int(v_['I+7']))])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.00001)
        self.xupper = np.full((self.n,1),+float('inf'))
        self.xupper[ix_['X'+str(int(v_['1']))]] = 2000.0
        self.xupper[ix_['X'+str(int(v_['2']))]] = 16000.0
        self.xupper[ix_['X'+str(int(v_['3']))]] = 120.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X'+str(int(v_['1']))]] = float(1745.0)
        self.x0[ix_['X'+str(int(v_['2']))]] = float(12000.0)
        self.x0[ix_['X'+str(int(v_['3']))]] = float(110.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eY2Y5', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY2', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY3', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY4', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY5', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY6', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY7', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        [it,iet_,_] = s2mpj_ii( 'eY8', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftv = loaset(elftv,it,2,'U3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'E25'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY2Y5')
        ielftype = arrset(ielftype,ie,iet_["eY2Y5"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY2')
        ielftype = arrset(ielftype,ie,iet_["eY2"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY3')
        ielftype = arrset(ielftype,ie,iet_["eY3"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY4')
        ielftype = arrset(ielftype,ie,iet_["eY4"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY5')
        ielftype = arrset(ielftype,ie,iet_["eY5"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY6')
        ielftype = arrset(ielftype,ie,iet_["eY6"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY7')
        ielftype = arrset(ielftype,ie,iet_["eY7"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eY8')
        ielftype = arrset(ielftype,ie,iet_["eY8"])
        vname = 'X1'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_] = s2mpj_nlx(self, vname,ix_,1,float(0.00001),None,None)
        posev = np.where(elftv[ielftype[ie]]=='U3')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E25'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(-0.063))
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        self.grelw = loaset(self.grelw,ig,posel,float(3.36))
        ig = ig_['AG'+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['1']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E2'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['2']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['3']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['4']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['4']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E5'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['5']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['5']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['6']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['6']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E7'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AG'+str(int(v_['7']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['AL'+str(int(v_['7']))]
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOI2-AN-3-14"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(self):

        import numpy as np
        self.efpar = np.array([]);
        self.efpar = arrset( self.efpar,0,.FALSE.0)
        return pbm

    @staticmethod
    def eY2Y5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if self.efpar[0]==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(2)*Y(5)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = Y(2)*G(5,1)+Y(5)*G(2,1)
            g_[1] = Y(2)*G(5,2)+Y(5)*G(2,2)
            g_[2] = Y(2)*G(5,3)+Y(5)*G(2,3)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = Y(2)*H(5,1)+2.0*G(5,1)*G(2,1)+Y(5)*H(2,1)
                H_[0,1] = Y(2)*H(5,2)+G(5,1)*G(2,2)+G(5,2)*G(2,1)+Y(5)*H(2,2)
                H_[1,0] = H_[0,1]
                H_[0,2] = Y(2)*H(5,3)+G(5,1)*G(2,3)+G(5,3)*G(2,1)+Y(5)*H(2,3)
                H_[2,0] = H_[0,2]
                H_[1,1] = Y(2)*H(5,4)+2.0*G(5,2)*G(2,2)+Y(5)*H(2,4)
                H_[1,2] = Y(2)*H(5,5)+G(5,2)*G(2,3)+G(5,3)*G(2,2)+Y(5)*H(2,5)
                H_[2,1] = H_[1,2]
                H_[2,2] = Y(2)*H(5,6)+2.0*G(5,3)*G(2,3)+Y(5)*H(2,6)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(2)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(2,1)
            g_[1] = G(2,2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = H(2,1)
                H_[1,0] = H(2,2)
                H_[0,1] = H_[1,0]
                H_[1,1] = H(2,4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(3)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(3,1)
            g_[1] = G(3,2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = H(3,1)
                H_[1,0] = H(3,2)
                H_[0,1] = H_[1,0]
                H_[1,1] = H(3,4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(4)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(4,1)
            g_[1] = G(4,2)
            g_[2] = G(4,3)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = H(4,1)
                H_[1,0] = H(4,2)
                H_[0,1] = H_[1,0]
                H_[2,0] = H(4,3)
                H_[0,2] = H_[2,0]
                H_[1,1] = H(4,4)
                H_[2,1] = H(4,5)
                H_[1,2] = H_[2,1]
                H_[2,2] = H(4,6)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY5(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(5)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(5,1)
            g_[1] = G(5,2)
            g_[2] = G(5,3)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = H(5,1)
                H_[1,0] = H(5,2)
                H_[0,1] = H_[1,0]
                H_[2,0] = H(5,3)
                H_[0,2] = H_[2,0]
                H_[1,1] = H(5,4)
                H_[2,1] = H(5,5)
                H_[1,2] = H_[2,1]
                H_[2,2] = H(5,6)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY6(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(6)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(6,1)
            g_[1] = G(6,2)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = H(6,1)
                H_[1,0] = H(6,2)
                H_[0,1] = H_[1,0]
                H_[1,1] = H(6,4)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY7(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        EVAL = .TRUE.0
        f_   = Y(7)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(7,1)
            g_[1] = G(7,2)
            g_[2] = G(7,3)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = H(7,1)
                H_[1,0] = H(7,2)
                H_[0,1] = H_[1,0]
                H_[2,0] = H(7,3)
                H_[0,2] = H_[2,0]
                H_[1,1] = H(7,4)
                H_[2,1] = H(7,5)
                H_[1,2] = H_[2,1]
                H_[2,2] = H(7,6)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY8(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        if EVAL==0:
            DUMMY = HS67(EV_[0],EV_[1],EV_[2],Y,G,H)
        f_   = Y(8)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G(8,1)
            g_[1] = G(8,2)
            g_[2] = G(8,3)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = H(8,1)
                H_[1,0] = H(8,2)
                H_[0,1] = H_[1,0]
                H_[2,0] = H(8,3)
                H_[0,2] = H_[2,0]
                H_[1,1] = H(8,4)
                H_[2,1] = H(8,5)
                H_[1,2] = H_[2,1]
                H_[2,2] = H(8,6)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

