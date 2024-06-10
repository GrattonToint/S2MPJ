from s2mpjlib import *
class  HANGING(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HANGING
#    *********
# 
#    A catenary problem in 3 dimensions.  A rectangular grid is hung from its
#    4 corners under gravity.  The problem is to determine the resulting shape.
# 
#    Source:  
#    an example in a talk by Nesterova and Vial, LLN, 1994.
# 
#    SIF input: Ph. Toint, November 1994.
# 
#    classification = "LQR2-AY-V-V"
# 
#    dimension of the grid
# 
#           Alternative values for the SIF file parameters:
# IE NX                  3              $-PARAMETER n = 27
# IE NY                  3              $-PARAMETER
# 
# IE NX                  5              $-PARAMETER n = 90
# IE NY                  6              $-PARAMETER
# 
# IE NX                  10             $-PARAMETER n = 300  original value
# IE NY                  10             $-PARAMETER
# 
# IE NX                  20             $-PARAMETER n = 1800
# IE NY                  30             $-PARAMETER
# 
# IE NX                  40             $-PARAMETER n = 3600
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HANGING'

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
            v_['NX'] = int(3);  #  SIF file default value
        else:
            v_['NX'] = int(args[0])
# IE NY                  30             $-PARAMETER
        if nargin<2:
            v_['NY'] = int(3);  #  SIF file default value
        else:
            v_['NY'] = int(args[1])
        v_['LX'] = 1.8
        v_['LY'] = 1.8
        v_['1'] = 1
        v_['NX-1'] = -1+v_['NX']
        v_['NY-1'] = -1+v_['NY']
        v_['LX2'] = v_['LX']*v_['LX']
        v_['LY2'] = v_['LY']*v_['LY']
        v_['RNX'] = float(v_['NX'])
        v_['RNY'] = float(v_['NY'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Y'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('Z'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Z'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['Z'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                [ig,ig_,_] = s2mpj_ii('RC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'RC'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('DC'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'DC'+str(I)+','+str(J))
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
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['RC'+str(I)+','+str(J)],float(v_['LX2']))
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['DC'+str(I)+','+str(J)],float(v_['LY2']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Z'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Z'+str(int(v_['1']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['NX']))+','+str(int(v_['1']))]] = v_['RNX']
        pb.xupper[ix_['X'+str(int(v_['NX']))+','+str(int(v_['1']))]] = v_['RNX']
        pb.xlower[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        pb.xupper[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['1']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        pb.xlower[ix_['Y'+str(int(v_['1']))+','+str(int(v_['NY']))]] = v_['RNY']
        pb.xupper[ix_['Y'+str(int(v_['1']))+','+str(int(v_['NY']))]] = v_['RNY']
        pb.xlower[ix_['Z'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        pb.xupper[ix_['Z'+str(int(v_['1']))+','+str(int(v_['NY']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNX']
        pb.xupper[ix_['X'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNX']
        pb.xlower[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNY']
        pb.xupper[ix_['Y'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = v_['RNY']
        pb.xlower[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = 0.0
        pb.xupper[ix_['Z'+str(int(v_['NX']))+','+str(int(v_['NY']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            v_['RI-1'] = float(v_['I-1'])
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                v_['RJ-1'] = float(v_['J-1'])
                pb.x0[ix_['X'+str(I)+','+str(J)]] = float(v_['RI-1'])
                pb.x0[ix_['Y'+str(I)+','+str(J)]] = float(v_['RJ-1'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftv = loaset(elftv,it,1,'YY')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for J in range(int(v_['1']),int(v_['NY-1'])+1):
            v_['J+1'] = 1+J
            for I in range(int(v_['1']),int(v_['NX'])+1):
                ename = 'RX'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'RY'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'RZ'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Z'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'DX'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'X'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'DY'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'Y'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'DZ'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                    ielftype = arrset( ielftype,ie,iet_['eISQ'])
                vname = 'Z'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XX')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Z'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='YY')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY-1'])+1):
                ig = ig_['RC'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RX'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RY'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RZ'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NX-1'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ig = ig_['DC'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DX'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DY'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['DZ'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(3,3)          -6.1184107487
# LO SOLTN(5,6)          -77.260229515
# LO SOLTN(10,10)        -620.17603242
# LO SOLTN(20,30)        -1025.4292887
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-AY-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

