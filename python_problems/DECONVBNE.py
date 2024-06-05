from s2mpjlib import *
class  DECONVBNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DECONVBNE
#    *********
# 
#    A problem arising in deconvolution analysis
#    (bounded variables version).
# 
#    Source:
#    J.P. Rasson, Private communication, 1996.
# 
#    SIF input: Ph. Toint, Nov 1996.
#    unititialized variables fixed at zero, Nick Gould, Feb, 2013
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "NOR2-MN-61-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DECONVBNE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DECONVBNE'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['LGSG'] = 11
        v_['LGTR'] = 40
        v_['-LGSG'] = -1*v_['LGSG']
        v_['PIC'] = 3.0000000000
        v_['TR1'] = 0.0000000000
        v_['TR2'] = 0.0000000000
        v_['TR3'] = 1.600000E-03
        v_['TR4'] = 5.400000E-03
        v_['TR5'] = 7.020000E-02
        v_['TR6'] = 0.1876000000
        v_['TR7'] = 0.3320000000
        v_['TR8'] = 0.7640000000
        v_['TR9'] = 0.9320000000
        v_['TR10'] = 0.8120000000
        v_['TR11'] = 0.3464000000
        v_['TR12'] = 0.2064000000
        v_['TR13'] = 8.300000E-02
        v_['TR14'] = 3.400000E-02
        v_['TR15'] = 6.179999E-02
        v_['TR16'] = 1.2000000000
        v_['TR17'] = 1.8000000000
        v_['TR18'] = 2.4000000000
        v_['TR19'] = 9.0000000000
        v_['TR20'] = 2.4000000000
        v_['TR21'] = 1.8010000000
        v_['TR22'] = 1.3250000000
        v_['TR23'] = 7.620000E-02
        v_['TR24'] = 0.2104000000
        v_['TR25'] = 0.2680000000
        v_['TR26'] = 0.5520000000
        v_['TR27'] = 0.9960000000
        v_['TR28'] = 0.3600000000
        v_['TR29'] = 0.2400000000
        v_['TR30'] = 0.1510000000
        v_['TR31'] = 2.480000E-02
        v_['TR32'] = 0.2432000000
        v_['TR33'] = 0.3602000000
        v_['TR34'] = 0.4800000000
        v_['TR35'] = 1.8000000000
        v_['TR36'] = 0.4800000000
        v_['TR37'] = 0.3600000000
        v_['TR38'] = 0.2640000000
        v_['TR39'] = 6.000000E-03
        v_['TR40'] = 6.000000E-03
        v_['SSG1'] = 1.000000E-02
        v_['SSG2'] = 2.000000E-02
        v_['SSG3'] = 0.4000000000
        v_['SSG4'] = 0.6000000000
        v_['SSG5'] = 0.8000000000
        v_['SSG6'] = 3.0000000000
        v_['SSG7'] = 0.8000000000
        v_['SSG8'] = 0.6000000000
        v_['SSG9'] = 0.4400000000
        v_['SSG10'] = 1.000000E-02
        v_['SSG11'] = 1.000000E-02
        v_['CC1'] = 0.0
        v_['CC2'] = 0.0
        v_['CC3'] = 0.0
        v_['CC4'] = 0.0
        v_['CC5'] = 0.0
        v_['CC6'] = 0.0
        v_['CC7'] = 0.0
        v_['CC8'] = 0.0
        v_['CC9'] = 0.0
        v_['CC10'] = 0.0
        v_['CC11'] = 0.0
        v_['CC12'] = 0.0
        v_['CC13'] = 0.0
        v_['CC14'] = 0.0
        v_['CC15'] = 0.0
        v_['CC16'] = 0.0
        v_['CC17'] = 0.0
        v_['CC18'] = 0.0
        v_['CC19'] = 0.0
        v_['CC20'] = 0.0
        v_['CC21'] = 0.0
        v_['CC22'] = 0.0
        v_['CC23'] = 0.0
        v_['CC24'] = 0.0
        v_['CC25'] = 0.0
        v_['CC26'] = 0.0
        v_['CC27'] = 0.0
        v_['CC28'] = 0.0
        v_['CC29'] = 0.0
        v_['CC30'] = 0.0
        v_['CC31'] = 0.0
        v_['CC32'] = 0.0
        v_['CC33'] = 0.0
        v_['CC34'] = 0.0
        v_['CC35'] = 0.0
        v_['CC36'] = 0.0
        v_['CC37'] = 0.0
        v_['CC38'] = 0.0
        v_['CC39'] = 0.0
        v_['CC40'] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for K in range(int(v_['-LGSG']),int(v_['LGTR'])+1):
            [iv,ix_,_] = s2mpj_ii('C'+str(K),ix_)
            pb.xnames=arrset(pb.xnames,iv,'C'+str(K))
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            [iv,ix_,_] = s2mpj_ii('SG'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'SG'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            [ig,ig_,_] = s2mpj_ii('R'+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(K))
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
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['R'+str(K)],float(v_['TR'+str(K)]))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            pb.xlower[ix_['SG'+str(I)]] = 0.0
            pb.xupper[ix_['SG'+str(I)]] = v_['PIC']
        for K in range(int(v_['-LGSG']),int(v_['0'])+1):
            pb.xlower[ix_['C'+str(K)]] = 0.0
            pb.xupper[ix_['C'+str(K)]] = 0.0
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            pb.xlower[ix_['C'+str(K)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            pb.x0[ix_['C'+str(K)]] = float(v_['CC'+str(K)])
        for I in range(int(v_['1']),int(v_['LGSG'])+1):
            pb.x0[ix_['SG'+str(I)]] = float(v_['SSG'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftp = []
        elftp = loaset(elftp,it,0,'IDX')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            for I in range(int(v_['1']),int(v_['LGSG'])+1):
                v_['K-I'] = K-I
                v_['K-I+1'] = 1+v_['K-I']
                v_['RIDX'] = float(v_['K-I+1'])
                ename = 'PROD'+str(K)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePR')
                ielftype = arrset(ielftype, ie, iet_["ePR"])
                vname = 'SG'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'C'+str(int(v_['K-I+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='IDX')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RIDX']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for K in range(int(v_['1']),int(v_['LGTR'])+1):
            for I in range(int(v_['1']),int(v_['LGSG'])+1):
                ig = ig_['R'+str(K)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PROD'+str(K)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-MN-61-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        NEGIDX = pbm.elpar[iel_][0]<=0.0
        if NEGIDX!=0:
            SCAL = 0.0
        if NEGIDX==0:
            SCAL = 1.0
        f_   = SCAL*EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = SCAL*EV_[1]
            g_[1] = SCAL*EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = SCAL
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

