from s2xlib import *
class  QR3D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QR3D
#    *********
# 
#    Find the QR factorization of a tridiagonal matrix A.
#    The problem is formulated as a system of quadratic equations
#    whose unknowns are the elements of the orthogonal matrix Q and of
#    the upper triangular matrix R.  In this version of the problem,
#    the banded structure of R is not imposed as a constraint. See problem
#    QR3DBD for the case where this structure is explicitly used.
# 
#    The problem is non-convex.
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, Nov 1993
# 
#    classification = "NQR2-AN-V-V"
# 
#    Define the matrix order M  ( M >= 3 ).
#    There are M * ( 3M + 1) / 2 variables and equations.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'QR3D'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'QR3D'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(5);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE M                   10             $-PARAMETER  n = 155  original value
# IE M                   20             $-PARAMETER  n = 610
        v_['1'] = 1
        v_['2'] = 2
        v_['M-1'] = -1+v_['M']
        v_['RM'] = float(v_['M'])
        v_['2/M'] = 2.0/v_['RM']
        v_['A'+str(int(v_['1']))+','+str(int(v_['1']))] = v_['2/M']
        v_['A'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.0
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['1-I'] = -1*v_['I-1']
            v_['R1-I'] = float(v_['1-I'])
            v_['1-I/M'] = v_['R1-I']/v_['RM']
            v_['2I'] = 2*I
            v_['R2I'] = float(v_['2I'])
            v_['2I/M'] = v_['R2I']/v_['RM']
            v_['A'+str(I)+','+str(int(v_['I-1']))] = v_['1-I/M']
            v_['A'+str(I)+','+str(I)] = v_['2I/M']
            v_['A'+str(I)+','+str(int(v_['I+1']))] = v_['1-I/M']
        v_['RM-1'] = float(v_['M-1'])
        v_['1-M'] = -1.0*v_['RM-1']
        v_['1-M/M'] = v_['1-M']/v_['RM']
        v_['2M'] = 2.0*v_['RM']
        v_['A'+str(int(v_['M']))+','+str(int(v_['M-1']))] = v_['1-M/M']
        v_['A'+str(int(v_['M']))+','+str(int(v_['M']))] = v_['2M']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [iv,ix_,_] = s2x_ii('Q'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'Q'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                [iv,ix_,_] = s2x_ii('R'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'R'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                [ig,ig_,_] = s2x_ii('O'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'O'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2x_ii('F'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'F'+str(I)+','+str(J))
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
            pbm.gconst = arrset(pbm.gconst,ig_['O'+str(I)+','+str(I)],float(1.0))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['F'+str(int(v_['1']))+','+str(int(v_['1']))],float(v_['A'+str(int(v_['1']))+','+str(int(v_['1']))])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['F'+str(int(v_['1']))+','+str(int(v_['2']))],float(v_['A'+str(int(v_['1']))+','+str(int(v_['2']))])))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['F'+str(I)+','+str(int(v_['I-1']))],float(v_['A'+str(I)+','+str(int(v_['I-1']))])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['F'+str(I)+','+str(I)],float(v_['A'+str(I)+','+str(I)])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['F'+str(I)+','+str(int(v_['I+1']))],float(v_['A'+str(I)+','+str(int(v_['I+1']))])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['F'+str(int(v_['M']))+','+str(int(v_['M-1']))],float(v_['A'+str(int(v_['M']))+','+str(int(v_['M-1']))])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['F'+str(int(v_['M']))+','+str(int(v_['M']))],float(v_['A'+str(int(v_['M']))+','+str(int(v_['M']))])))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pb.xlower[ix_['R'+str(I)+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            pb.x0[ix_['Q'+str(I)+','+str(I)]] = float(1.0)
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            v_['I+1'] = 1+I
            pb.x0[ix_['R'+str(I)+','+str(I)]] = float(v_['A'+str(I)+','+str(I)])
            pb.x0[ix_['R'+str(I)+','+str(int(v_['I+1']))]]  = (
                  float(v_['A'+str(I)+','+str(int(v_['I+1']))]))
        pb.x0[ix_['R'+str(int(v_['M']))+','+str(int(v_['M']))]]  = (
              float(v_['A'+str(int(v_['M']))+','+str(int(v_['M']))]))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    ename = 'C'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2x_ii(ename,ie_)
                    if newelt:
                        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                        ielftype = arrset( ielftype,ie,iet_['en2PR'])
                    vname = 'Q'+str(I)+','+str(K)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'Q'+str(J)+','+str(K)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for K in range(int(v_['1']),int(J)+1):
                    ename = 'B'+str(I)+','+str(J)+','+str(K)
                    [ie,ie_,newelt] = s2x_ii(ename,ie_)
                    if newelt:
                        pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                        ielftype = arrset( ielftype,ie,iet_['en2PR'])
                    vname = 'Q'+str(I)+','+str(K)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'R'+str(K)+','+str(J)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(I),int(v_['M'])+1):
                for K in range(int(v_['1']),int(v_['M'])+1):
                    ig = ig_['O'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['C'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for K in range(int(v_['1']),int(J)+1):
                    ig = ig_['F'+str(I)+','+str(J)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B'+str(I)+','+str(J)+','+str(K)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NQR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

