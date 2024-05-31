from s2xlib import *
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
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TWOD'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'TWOD'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(2);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE N                   40             $-PARAMETER
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                for K in range(int(v_['0']),int(v_['M'])+1):
                    [iv,ix_,_] = s2x_ii('Y'+str(K)+','+str(I)+','+str(J),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'Y'+str(K)+','+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2x_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N1'])+1):
            v_['I+'] = 1+I
            v_['I-'] = -1+I
            for J in range(int(v_['1']),int(v_['N1'])+1):
                v_['J+'] = 1+J
                v_['J-'] = -1+J
                for K in range(int(v_['0']),int(v_['M1'])+1):
                    v_['K+'] = 1+K
                    [ig,ig_,_] = s2x_ii('P'+str(K)+','+str(I)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'P'+str(K)+','+str(I)+','+str(J))
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(J)]
                    pbm.A[ig,iv] = float(v_['1/DT+2/H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(J)]
                    pbm.A[ig,iv] = float(v_['-1/DT+2/H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['J-']))]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['J+']))]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(int(v_['I-']))+','+str(J)]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(K)+','+str(int(v_['I+']))+','+str(J)]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(int(v_['I-']))+','+str(J)]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(int(v_['I+']))+','+str(J)]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(int(v_['J-']))]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
                    iv = ix_['Y'+str(int(v_['K+']))+','+str(I)+','+str(int(v_['J+']))]
                    pbm.A[ig,iv] = float(v_['-.1/2H2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N1'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2x_ii('B1'+str(K)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B1'+str(K)+','+str(I))
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N2']))]
                pbm.A[ig,iv] = float(v_['1/2DY'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N1']))]
                pbm.A[ig,iv] = float(v_['-2/DY'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['N']))]
                pbm.A[ig,iv] = float(v_['3/2DY+1'])+pbm.A[ig,iv]
                iv = ix_['U'+str(K)+','+str(I)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('B2'+str(K)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B2'+str(K)+','+str(I))
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['2']))]
                pbm.A[ig,iv] = float(v_['1/2DY'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['1']))]
                pbm.A[ig,iv] = float(v_['-2/DY'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(I)+','+str(int(v_['0']))]
                pbm.A[ig,iv] = float(v_['3/2DY+1'])+pbm.A[ig,iv]
        for J in range(int(v_['1']),int(v_['N1'])+1):
            for K in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2x_ii('B3'+str(K)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B3'+str(K)+','+str(J))
                iv = ix_['Y'+str(K)+','+str(int(v_['2']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/2DX'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/DX'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['0']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['3/2DX+1'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2x_ii('B4'+str(K)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B4'+str(K)+','+str(J))
                iv = ix_['Y'+str(K)+','+str(int(v_['N2']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/2DX'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['N1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/DX'])+pbm.A[ig,iv]
                iv = ix_['Y'+str(K)+','+str(int(v_['N']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['3/2DX+1'])+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                pb.xlower[ix_['Y'+str(int(v_['0']))+','+str(I)+','+str(J)]] = 0.0
                pb.xupper[ix_['Y'+str(int(v_['0']))+','+str(I)+','+str(J)]] = 0.0
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    pb.xlower[ix_['Y'+str(K)+','+str(I)+','+str(J)]] = 0.0
                    pb.xupper[ix_['Y'+str(K)+','+str(I)+','+str(J)]] = 0.8
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                pb.xlower[ix_['U'+str(I)+','+str(J)]] = 0.0
                pb.xupper[ix_['U'+str(I)+','+str(J)]] = v_['UA']
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        pb.y0 = np.full((pb.m,1),float(0.0))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                if('U'+str(I)+','+str(J) in ix_):
                    pb.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['UA'])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['U'+str(I)+','+str(J)]),float(v_['UA'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        [it,iet_,_] = s2x_ii( 'eSQD', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'YP')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['.5DXDYI'] = v_['.5DXDY']*v_['RI']
            for J in range(int(v_['0']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['.5DXDYIJ'] = v_['.5DXDYI']*v_['RJ']
                v_['YP'] = 0.25+v_['.5DXDYIJ']
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQD')
                ielftype = arrset(ielftype, ie, iet_["eSQD"])
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                vname = 'Y'+str(int(v_['M']))+','+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'E'+str(int(v_['M']))+','+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                posep = find(elftp[ielftype[ie]],lambda x:x=='YP')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['YP']))
        for K in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N1'])+1):
                ename = 'E'+str(K)+','+str(I)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'U'+str(K)+','+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(int(v_['N']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['0']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.125DXDY']))
        posel = len(pbm.grelt[ig])
        pbm.grelt  = (
              loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['N']))]))
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.125DXDY']))
        for J in range(int(v_['1']),int(v_['N1'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['0']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(J)]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(J)+','+str(int(v_['0']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.25DXDY']))
            posel = len(pbm.grelt[ig])
            pbm.grelt  = (
                  loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(int(v_['N']))+','+str(int(v_['N']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.25DXDY']))
        for I in range(int(v_['1']),int(v_['N1'])+1):
            for J in range(int(v_['1']),int(v_['N1'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(I)+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.5DXDY']))
        for K in range(int(v_['1']),int(v_['M1'])+1):
            for I in range(int(v_['1']),int(v_['N1'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(K)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.5ADTDX']))
        for I in range(int(v_['1']),int(v_['N1'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(int(v_['M']))+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['.25ADTDX']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(pbm,nargout,*args):

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
    def eSQD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = (EV_[0]-pbm.elpar[iel_][0])**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 2.0*(EV_[0]-pbm.elpar[iel_][0])
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

