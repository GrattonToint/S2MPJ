from s2xlib import *
class  HS86(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS86
#    *********
# 
#    Source: problem 86 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "OLR2-AN-5-10"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS86'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS86'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 5
        v_['M'] = 10
        v_['1'] = 1
        v_['E1'] = -15.0
        v_['E2'] = -27.0
        v_['E3'] = -36.0
        v_['E4'] = -18.0
        v_['E5'] = -12.0
        v_['C1,1'] = 30.0
        v_['C2,1'] = -20.0
        v_['C3,1'] = -10.0
        v_['C4,1'] = 32.0
        v_['C5,1'] = -10.0
        v_['C1,2'] = -20.0
        v_['C2,2'] = 39.0
        v_['C3,2'] = -6.0
        v_['C4,2'] = -31.0
        v_['C5,2'] = 32.0
        v_['C1,3'] = -10.0
        v_['C2,3'] = -6.0
        v_['C3,3'] = 10.0
        v_['C4,3'] = -6.0
        v_['C5,3'] = -10.0
        v_['C1,4'] = 32.0
        v_['C2,4'] = -31.0
        v_['C3,4'] = -6.0
        v_['C4,4'] = 39.0
        v_['C5,4'] = -20.0
        v_['C1,5'] = -10.0
        v_['C2,5'] = 32.0
        v_['C3,5'] = -10.0
        v_['C4,5'] = -20.0
        v_['C5,5'] = 30.0
        v_['D1'] = 4.0
        v_['D2'] = 8.0
        v_['D3'] = 10.0
        v_['D4'] = 6.0
        v_['D5'] = 2.0
        v_['A1,1'] = -16.0
        v_['A2,1'] = 0.0
        v_['A3,1'] = -3.5
        v_['A4,1'] = 0.0
        v_['A5,1'] = 0.0
        v_['A6,1'] = 2.0
        v_['A7,1'] = -1.0
        v_['A8,1'] = -1.0
        v_['A9,1'] = 1.0
        v_['A10,1'] = 1.0
        v_['A1,2'] = 2.0
        v_['A2,2'] = -2.0
        v_['A3,2'] = 0.0
        v_['A4,2'] = -2.0
        v_['A5,2'] = -9.0
        v_['A6,2'] = 0.0
        v_['A7,2'] = -1.0
        v_['A8,2'] = -2.0
        v_['A9,2'] = 2.0
        v_['A10,2'] = 1.0
        v_['A1,3'] = 0.0
        v_['A2,3'] = 0.0
        v_['A3,3'] = 2.0
        v_['A4,3'] = 0.0
        v_['A5,3'] = -2.0
        v_['A6,3'] = -4.0
        v_['A7,3'] = -1.0
        v_['A8,3'] = -3.0
        v_['A9,3'] = 3.0
        v_['A10,3'] = 1.0
        v_['A1,4'] = 1.0
        v_['A2,4'] = 4.0
        v_['A3,4'] = 0.0
        v_['A4,4'] = -4.0
        v_['A5,4'] = 1.0
        v_['A6,4'] = 0.0
        v_['A7,4'] = -1.0
        v_['A8,4'] = -2.0
        v_['A9,4'] = 4.0
        v_['A10,4'] = 1.0
        v_['A1,5'] = 0.0
        v_['A2,5'] = 2.0
        v_['A3,5'] = 0.0
        v_['A4,5'] = -1.0
        v_['A5,5'] = -2.8
        v_['A6,5'] = 0.0
        v_['A7,5'] = -1.0
        v_['A8,5'] = -1.0
        v_['A9,5'] = 5.0
        v_['A10,5'] = 1.0
        v_['B1'] = -40.0
        v_['B2'] = -2.0
        v_['B3'] = -0.25
        v_['B4'] = -4.0
        v_['B5'] = -4.0
        v_['B6'] = -1.0
        v_['B7'] = -40.0
        v_['B8'] = -60.0
        v_['B9'] = 5.0
        v_['B10'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2x_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(J)]
            pbm.A[ig,iv] = float(v_['E'+str(J)])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2x_ii('C'+str(I),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(I))
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['A'+str(I)+','+str(J)])+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)],float(v_['B'+str(I)]))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        pb.y0 = np.full((pb.m,1),float(0.0))
        if('X5' in ix_):
            pb.x0[ix_['X5']] = float(1.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X5']),float(1.0)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eCUBE', iet_)
        elftv = loaset(elftv,it,0,'XJ')
        [it,iet_,_] = s2x_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'XI')
        elftv = loaset(elftv,it,1,'XJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for J in range(int(v_['1']),int(v_['N'])+1):
            ename = 'D'+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eCUBE')
            ielftype = arrset(ielftype, ie, iet_["eCUBE"])
            vname = 'X'+str(J)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='XJ')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            for I in range(int(v_['1']),int(v_['N'])+1):
                ename = 'C'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XJ')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['D'+str(J)]))
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['C'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
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
        pb.pbclass = "OLR2-AN-5-10"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eCUBE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0e+0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0e+0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(pbm,nargout,*args):

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
                H_[0,1] = 1.0e+0
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
