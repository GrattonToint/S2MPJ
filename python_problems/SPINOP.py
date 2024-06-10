from s2mpjlib import *
class  SPINOP(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPINOP
#    *********
# 
#    Given n particles z_j = x_j + i * y_j in the complex plane, 
#    determine their positions so that the equations
# 
#      z'_j = lambda z_j, 
# 
#    where z_j = sum_k \j i / conj( z_j - z_k ) and i = sqrt(-1)
#    for some lamda = mu + i * omega
# 
#    Pick the solution for which sum_i sum_k\i |z_i -z_j|^2 is smllest
# 
#    A problem posed by Nick Trefethen
# 
#    SIF input: Nick Gould, June 2009
# 
#    classification = "QOR2-AN-V-V"
# 
#    Number of particles n
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER matrix dimension
# IE N                   3              $-PARAMETER matrix dimension
# IE N                   5              $-PARAMETER matrix dimension
# IE N                   10             $-PARAMETER matrix dimension
# IE N                   15             $-PARAMETER matrix dimension
# IE N                   20             $-PARAMETER matrix dimension
# IE N                   25             $-PARAMETER matrix dimension
# IE N                   30             $-PARAMETER matrix dimension
# IE N                   35             $-PARAMETER matrix dimension
# IE N                   40             $-PARAMETER matrix dimension
# IE N                   45             $-PARAMETER matrix dimension
# IE N                   50             $-PARAMETER matrix dimension
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPINOP'

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
            v_['N'] = int(2);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['1'] = 1
        v_['2'] = 2
        v_['RN'] = float(v_['N'])
        v_['PI/4'] = np.arctan(1.0)
        v_['2PI'] = 8.0*v_['PI/4']
        v_['2PI/N'] = v_['2PI']/v_['RN']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('MU',ix_)
        pb.xnames=arrset(pb.xnames,iv,'MU')
        [iv,ix_,_] = s2mpj_ii('OMEGA',ix_)
        pb.xnames=arrset(pb.xnames,iv,'OMEGA')
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('R'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('I'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'I'+str(I))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                [ig,ig_,_] = s2mpj_ii('M'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'M'+str(I)+','+str(J))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['2PII/N'] = v_['2PI/N']*v_['RI']
            v_['C'] = np.cos(v_['2PII/N'])
            v_['S'] = np.sin(v_['2PII/N'])
            pb.x0[ix_['X'+str(I)]] = float(v_['C'])
            pb.x0[ix_['Y'+str(I)]] = float(v_['S'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eRATIO', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        elftv = loaset(elftv,it,2,'V')
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eSQR2', iet_)
        elftv = loaset(elftv,it,0,'Y')
        elftv = loaset(elftv,it,1,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'MX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'MU'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'MY'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'MU'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'OX'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'OMEGA'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'OY'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
            ielftype = arrset(ielftype, ie, iet_["en2PR"])
            vname = 'Y'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'OMEGA'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ename = 'RX'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eRATIO')
                ielftype = arrset(ielftype, ie, iet_["eRATIO"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'RY'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eRATIO')
                ielftype = arrset(ielftype, ie, iet_["eRATIO"])
                vname = 'Y'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ename = 'V'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
                ielftype = arrset(ielftype, ie, iet_["eSQR"])
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'X'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQR2')
                ielftype = arrset(ielftype, ie, iet_["eSQR2"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Y'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQR2')
                ielftype = arrset(ielftype, ie, iet_["eSQR2"])
                vname = 'Y'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Y'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['V'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            ig = ig_['R'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['MX'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OY'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['R'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RY'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['R'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RY'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            ig = ig_['I'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['MY'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['OX'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['I'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RX'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['I'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RX'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['M'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['V'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['X'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQR(pbm,nargout,*args):

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

    @staticmethod
    def eSQR2(pbm,nargout,*args):

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

    @staticmethod
    def en2PR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
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

    @staticmethod
    def eRATIO(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,2] = U_[1,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]/IV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.0/IV_[1]**2
            g_[1] = -2.0*IV_[0]/IV_[1]**3
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -2.0/IV_[1]**3
                H_[1,0] = H_[0,1]
                H_[1,1] = 6.0*IV_[0]/IV_[1]**4
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

