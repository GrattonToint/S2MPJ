from s2mpjlib import *
class  HS84(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    *********
#    Problem : HS84
#    *********
# 
#    Source: problem 84 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "C-CQQR2-AN-5-3"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS84'

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
        v_['15'] = 15
        v_['16'] = 16
        v_['17'] = 17
        v_['18'] = 18
        v_['19'] = 19
        v_['20'] = 20
        v_['21'] = 21
        v_['A'+str(int(v_['1']))] = -24345.0
        v_['A'+str(int(v_['2']))] = -8720288.849
        v_['MA'+str(int(v_['2']))] = -1.0*v_['A'+str(int(v_['2']))]
        v_['A'+str(int(v_['3']))] = 150512.5253
        v_['MA'+str(int(v_['3']))] = -1.0*v_['A'+str(int(v_['3']))]
        v_['A'+str(int(v_['4']))] = -156.6950325
        v_['MA'+str(int(v_['4']))] = -1.0*v_['A'+str(int(v_['4']))]
        v_['A'+str(int(v_['5']))] = 476470.3222
        v_['MA'+str(int(v_['5']))] = -1.0*v_['A'+str(int(v_['5']))]
        v_['A'+str(int(v_['6']))] = 729482.8271
        v_['MA'+str(int(v_['6']))] = -1.0*v_['A'+str(int(v_['6']))]
        v_['A'+str(int(v_['7']))] = -145421.402
        v_['A'+str(int(v_['8']))] = 2931.1506
        v_['A'+str(int(v_['9']))] = -40.427932
        v_['A'+str(int(v_['10']))] = 5106.192
        v_['A'+str(int(v_['11']))] = 15711.36
        v_['A'+str(int(v_['12']))] = -155011.1084
        v_['A'+str(int(v_['13']))] = 4360.53352
        v_['A'+str(int(v_['14']))] = 12.9492344
        v_['A'+str(int(v_['15']))] = 10236.884
        v_['A'+str(int(v_['16']))] = 13176.786
        v_['A'+str(int(v_['17']))] = -326669.5104
        v_['A'+str(int(v_['18']))] = 7390.68412
        v_['A'+str(int(v_['19']))] = -27.8986976
        v_['A'+str(int(v_['20']))] = 16643.076
        v_['A'+str(int(v_['21']))] = 30988.146
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['5'])+1):
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
        valA = np.append(valA,float(v_['MA'+str(int(v_['2']))]))
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON1')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(v_['A'+str(int(v_['7']))]))
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON2')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(v_['A'+str(int(v_['12']))]))
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'CON3')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['X'+str(int(v_['1']))]])
        valA = np.append(valA,float(v_['A'+str(int(v_['17']))]))
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
        self.gconst  = (
              arrset(self.gconst,ig_['OBJ'],float(v_['A'+str(int(v_['1']))])))
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = np.full((ngrp,1),None)
        grange[gegrps] = np.full((self.nge,1),float('inf'))
        grange = arrset(grange,ig_['CON1'],float(294000.0))
        grange = arrset(grange,ig_['CON2'],float(294000.0))
        grange = arrset(grange,ig_['CON3'],float(277200.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X'+str(int(v_['1']))]] = 0.0
        self.xlower[ix_['X'+str(int(v_['2']))]] = 1.2
        self.xlower[ix_['X'+str(int(v_['3']))]] = 20.0
        self.xlower[ix_['X'+str(int(v_['4']))]] = 9.0
        self.xlower[ix_['X'+str(int(v_['5']))]] = 6.5
        self.xupper[ix_['X'+str(int(v_['1']))]] = 1000.0
        self.xupper[ix_['X'+str(int(v_['2']))]] = 2.4
        self.xupper[ix_['X'+str(int(v_['3']))]] = 60.0
        self.xupper[ix_['X'+str(int(v_['4']))]] = 9.3
        self.xupper[ix_['X'+str(int(v_['5']))]] = 7.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X'+str(int(v_['1']))]] = float(2.52)
        self.x0[ix_['X'+str(int(v_['2']))]] = float(2.0)
        self.x0[ix_['X'+str(int(v_['3']))]] = float(37.5)
        self.x0[ix_['X'+str(int(v_['4']))]] = float(9.25)
        self.x0[ix_['X'+str(int(v_['5']))]] = float(6.8)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['4'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            v_['IP1'] = 1+I
            vname = 'X'+str(int(v_['1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='U1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['IP1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='U2')[0]
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
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['MA'+str(int(v_['3']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['MA'+str(int(v_['4']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['MA'+str(int(v_['5']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['MA'+str(int(v_['6']))]))
        ig = ig_['CON1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['8']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['9']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['10']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['11']))]))
        ig = ig_['CON2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['13']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['14']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['15']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['16']))]))
        ig = ig_['CON3']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['1']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['18']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['2']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['19']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['3']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['20']))]))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['4']))])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['A'+str(int(v_['21']))]))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        self.cupper[np.arange(self.nle+self.neq,self.m)] = grange[gegrps]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQQR2-AN-5-3"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0,0]*EV_[1,0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]
            g_[1] = EV_[0,0]
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

