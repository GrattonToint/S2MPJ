from s2mpjlib import *
class  TWOD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TWOD
#    *********
# 
#    The twod_0 & _00.mod AMPL models from Hans Mittelmann (mittelmann@asu.edu)
#    See: http://plato.asu.edu/ftp/barrier/
# 
#    SIF input: Nick Gould, April 25th 2012
# 
#    classification = "QLR2-AN-V-V"
# 
#    the x-y discretization 
# 
#           Alternative values for the SIF file parameters:
# IE N                    2             $-PARAMETER
# IE N                   40             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TWOD'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(2);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   79             $-PARAMETER     twod_000.mod value
# IE N                   99             $-PARAMETER     twod_0.mod value
        v_['M'] = v_['N']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['ONE'] = 1.0
        v_['HALF'] = 0.5
        v_['-HALF'] = -0.5
        v_['A'] = 0.001
        v_['UA'] = 2.0
        v_['N1'] = -1+v_['N']
        v_['N2'] = -2+v_['N']
        v_['M1'] = -1+v_['M']
        v_['RN'] = float(v_['N'])
        v_['RM'] = float(v_['M'])
        v_['DX'] = v_['ONE']/v_['RN']
        v_['DY'] = v_['ONE']/v_['RM']
        v_['T'] = v_['ONE']
        v_['DT'] = v_['T']/v_['RM']
        v_['H2'] = v_['DX']*v_['DX']
        v_['DXDY'] = v_['DX']*v_['DY']
        v_['.5DXDY'] = 0.5*v_['DXDY']
        v_['.25DXDY'] = 0.25*v_['DXDY']
        v_['.125DXDY'] = 0.125*v_['DXDY']
        v_['DTDX'] = v_['DT']*v_['DX']
        v_['ADTDX'] = v_['A']*v_['DTDX']
        v_['.5ADTDX'] = 0.5*v_['ADTDX']
        v_['.25ADTDX'] = 0.5*v_['ADTDX']
        v_['1/2DX'] = v_['HALF']/v_['DX']
        v_['3/2DX'] = 3.0*v_['1/2DX']
        v_['-2/DX'] = -4.0*v_['1/2DX']
        v_['3/2DX+1'] = 1.0+v_['3/2DX']
        v_['1/2DY'] = v_['HALF']/v_['DY']
        v_['3/2DY'] = 3.0*v_['1/2DY']
        v_['-2/DY'] = -4.0*v_['1/2DY']
        v_['3/2DY+1'] = 1.0+v_['3/2DY']
        v_['1/DT'] = v_['ONE']/v_['DT']
        v_['-1/DT'] = -1.0*v_['1/DT']
        v_['-.1/2H2'] = v_['-HALF']/v_['H2']
        v_['2/H2'] = -4.0*v_['-.1/2H2']
        v_['1/DT+2/H2'] = v_['1/DT']+v_['2/H2']
        v_['-1/DT+2/H2'] = v_['-1/DT']+v_['2/H2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                for K in range(int(v_['0']),int(v_['M'])+1):
                    [iv,ix_,_] = s2mpj_ii('Y'+str(K)+','+str(I)+','+str(J),ix_)
                    self.xnames=arrset(self.xnames,iv,'Y'+str(K)+','+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N1'])+1):
            v_['I+'] = 1+I
            v_['I-'] = -1+I
            for J in range(int(v_['1']),int(v_['N1'])+1):
                v_['J+'] = 1+J
                v_['J-'] = -1+J
                for K in range(int(v_['0']),int(v_['M1'])+1):
                    v_['K+'] = 1+K
                    [ig,ig_,_] = s2mpj_ii('P'+str(K)+','+str(I)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'P'+str(K)+','+str(I)+','+str(J))
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(J)]
                    self.A[ig,iv] = float(v_['1/DT+2/H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(J)]
                    self.A[ig,iv] = float(v_['-1/DT+2/H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['J-']))]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['J+']))]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(int(v_['I-']))+','+str(J)]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(int(v_['I+']))+','+str(J)]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(int(v_['I-']))+','+str(J)]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(int(v_['I+']))+','+str(J)]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(int(v_['J-']))]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(int(v_['J+']))]
                    self.A[ig,iv] = float(v_['-.1/2H2'])+self.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N1'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('B1'+str(K)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B1'+str(K)+','+str(I))
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N2']))]
                self.A[ig,iv] = float(v_['1/2DY'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N1']))]
                self.A[ig,iv] = float(v_['-2/DY'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N']))]
                self.A[ig,iv] = float(v_['3/2DY+1'])+self.A[ig,iv]
                iv = ix_['U'+str(K)+','+str(I)]
                self.A[ig,iv] = float(-1.0)+self.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('B2'+str(K)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B2'+str(K)+','+str(I))
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['2']))]
                self.A[ig,iv] = float(v_['1/2DY'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['1']))]
                self.A[ig,iv] = float(v_['-2/DY'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['0']))]
                self.A[ig,iv] = float(v_['3/2DY+1'])+self.A[ig,iv]
        for J in range(int(v_['1']),int(v_['N1'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('B3'+str(K)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B3'+str(K)+','+str(J))
                iv = ix_['Y'+str(K)+','+str(int(v_['2']))+','+str(J)]
                self.A[ig,iv] = float(v_['1/2DX'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['1']))+','+str(J)]
                self.A[ig,iv] = float(v_['-2/DX'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['0']))+','+str(J)]
                self.A[ig,iv] = float(v_['3/2DX+1'])+self.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('B4'+str(K)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B4'+str(K)+','+str(J))
                iv = ix_['Y'+str(K)+','+str(int(v_['N2']))+','+str(J)]
                self.A[ig,iv] = float(v_['1/2DX'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['N1']))+','+str(J)]
                self.A[ig,iv] = float(v_['-2/DX'])+self.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['N']))+','+str(J)]
                self.A[ig,iv] = float(v_['3/2DX+1'])+self.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                self.xlower[ix_['Y'+str(int(v_['0']))+','+str(I)+','+str(J)]] = 0.0
                self.xupper[ix_['Y'+str(int(v_['0']))+','+str(I)+','+str(J)]] = 0.0
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    self.xlower[ix_['Y'+str(K)+','+str(I)+','+str(J)]] = 0.0
                    self.xupper[ix_['Y'+str(K)+','+str(I)+','+str(J)]] = 0.8
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                self.xlower[ix_['U'+str(I)+','+str(J)]] = 0.0
                self.xupper[ix_['U'+str(I)+','+str(J)]] = v_['UA']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        self.y0 = np.full((self.m,1),float(0.0))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                if('U'+str(I)+','+str(J) in ix_):
                    self.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['UA'])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)]),float(v_['UA'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        [it,iet_,_] = s2mpj_ii( 'eSQD', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'YP')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['.5DXDYI'] = v_['.5DXDY']*v_['RI']
            for J in range(int(v_['0']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['.5DXDYIJ'] = v_['.5DXDYI']*v_['RJ']
                v_['YP'] = 0.25+v_['.5DXDYIJ']
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQD')
                ielftype = arrset(ielftype, ie, iet_["eSQD"])
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                vname = 'Y'+str(int(v_['M']))+','+str(I)+','+str(J)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,0.0)
                posev = np.where(elftv[ielftype[ie]]=='Y')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                posep = np.where(elftp[ielftype[ie]]=='YP')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['YP']))
        for K in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N1'])+1):
                ename = 'E'+str(K)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'U'+str(K)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,0.0)
                posev = np.where(elftv[ielftype[ie]]=='U')[0]
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
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(int(v_['N']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(self.grelt[ig])
        self.grelt  = (
              loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['N']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(v_['.125DXDY']))
        for J in range(int(v_['1']),int(v_['N1'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(J)+','+str(int(v_['0']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['N']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['.25DXDY']))
        for I in range(int(v_['1']),int(v_['N1'])+1):
            for J in range(int(v_['1']),int(v_['N1'])+1):
                ig = ig_['OBJ']
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(I)+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['.5DXDY']))
        for K in range(int(v_['1']),int(v_['M1'])+1):
            for I in range(int(v_['1']),int(v_['N1'])+1):
                ig = ig_['OBJ']
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(K)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['.5ADTDX']))
        for I in range(int(v_['1']),int(v_['N1'])+1):
            ig = ig_['OBJ']
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['.25ADTDX']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        self.A.resize(ngrp,self.n)
        self.A     = self.A.tocsr()
        sA1,sA2    = self.A.shape
        self.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons =  np.where(self.congrps in np.setdiff1d(nlc,self.congrps))[0]
        self.pbclass = "QLR2-AN-V-V"

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
            g_[0] = 2.0*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSQD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-self.elpar[iel_][0])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-self.elpar[iel_][0])
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

