from s2mpjlib import *
class  LIPPERT2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LIPPERT2
#    *********
# 
#    A discrete approximation to a continuum optimal flow problem
#    in the unit square. The continuum problem requires that the
#    divergence of a given flow should be given everywhere in the
#    region of interest, with the restriction that the capacity of
#    the flow is bounded. The aim is then to maximize the given flow.
# 
#    The discrete problem (dual formulation 2) in the unit square is to 
#      minimize r
#      subject to dx( u_ij - ui-1j ) + dx( v_ij - vij-1 ) = s_ij
#                 u_ij^2 + v_ij^2 <= r^2
#                 u_i-1j^2 + v_ij^2 <= r^2
#                 u_ij^2 + v_ij-1^2 <= r^2
#                 u_i-1j^2 + v_ij-1^2 <= r^2
#      where 1 <= i <= nx, 1 <= j <= ny
#      and        r >= 0
# 
#    Source: R. A. Lippert
#      "Discrete approximations to continuum optimal flow problems"
#      Tech. Report, Dept of Maths, M.I.T., 2006
#    following a suggestion by Gil Strang
# 
#    SIF input: Nick Gould, September 2006
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "LQR2-MN-V-V"
# 
#    Number of nodes in x direction
# 
#           Alternative values for the SIF file parameters:
# IE NX                  2              $-PARAMETER
# IE NX                  3              $-PARAMETER
# IE NX                  10             $-PARAMETER
# IE NX                  40             $-PARAMETER
# IE NX                  100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LIPPERT2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'LIPPERT2'
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
# IE NY                  2              $-PARAMETER
# IE NY                  3              $-PARAMETER
# IE NY                  10             $-PARAMETER 
# IE NY                  40             $-PARAMETER
# IE NY                  100            $-PARAMETER
        if nargin<2:
            v_['NY'] = int(10);  #  SIF file default value
        else:
            v_['NY'] = int(args[1])
        v_['X+'] = 1+v_['NX']
        v_['X-'] = -1+v_['NX']
        v_['Y+'] = 1+v_['NY']
        v_['Y-'] = -1+v_['NY']
        v_['1'] = 1
        v_['0'] = 0
        v_['HALF'] = 0.5
        v_['ONE'] = 1.0
        v_['-ONE'] = -1.0
        v_['S'] = 1.0
        v_['-S'] = v_['S']*v_['-ONE']
        v_['RX'] = float(v_['NX'])
        v_['DX'] = v_['ONE']/v_['RX']
        v_['-DX'] = v_['-ONE']/v_['RX']
        v_['DX/2'] = v_['DX']*v_['HALF']
        v_['RY'] = float(v_['NY'])
        v_['DY'] = v_['ONE']/v_['RY']
        v_['-DY'] = v_['-ONE']/v_['RY']
        v_['DY/2'] = v_['DY']*v_['HALF']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('R',ix_)
        pb.xnames=arrset(pb.xnames,iv,'R')
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
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
        iv = ix_['R']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                [ig,ig_,_] = s2mpj_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'O'+str(I)+','+str(J))
                iv = ix_['U'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['DX'])+pbm.A[ig,iv]
                iv = ix_['U'+str(int(v_['I-1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['-DX'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['DY'])+pbm.A[ig,iv]
                iv = ix_['V'+str(I)+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(v_['-DY'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'B'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'C'+str(I)+','+str(J))
                [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'D'+str(I)+','+str(J))
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
            for J in range(int(v_['1']),int(v_['NY'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)+','+str(J)],float(v_['S']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower[ix_['R']] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        pb.x0[ix_['R']] = float(1.0)
        v_['ALPHA'] = 0.0
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                pb.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['ALPHA'])
            v_['ALPHA'] = v_['ALPHA']+v_['DX/2']
        v_['ALPHA'] = 0.0
        for J in range(int(v_['0']),int(v_['NY'])+1):
            for I in range(int(v_['1']),int(v_['NX'])+1):
                pb.x0[ix_['V'+str(I)+','+str(J)]] = float(v_['ALPHA'])
            v_['ALPHA'] = v_['ALPHA']+v_['DX/2']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQR', iet_)
        elftv = loaset(elftv,it,0,'ALPHA')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'RHO2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
            ielftype = arrset( ielftype,ie,iet_['eSQR'])
        vname = 'R'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='ALPHA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['1']),int(v_['NY'])+1):
                ename = 'P'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
                    ielftype = arrset( ielftype,ie,iet_['eSQR'])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='ALPHA')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                ename = 'Q'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eSQR')
                    ielftype = arrset( ielftype,ie,iet_['eSQR'])
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='ALPHA')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NX'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['NY'])+1):
                v_['J-1'] = -1+J
                ig = ig_['A'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['B'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['P'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['C'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                ig = ig_['D'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['P'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['Q'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['RHO2'])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               3.77245385
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LQR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

