from s2xlib import *
class  BRATU2D(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BRATU2D
#    *********
# 
#    The 2D Bratu problem on the unit square, using finite differences.
# 
#    Source: problem 3 in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "NOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit square.
#    There are P*P variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BRATU2D'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'BRATU2D'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(7);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
#           Alternative values for the SIF file parameters:
# IE P                   10             $-PARAMETER  n=P**2
# IE P                   22             $-PARAMETER  n=P**2
# IE P                   32             $-PARAMETER  n=P**2
# IE P                   72             $-PARAMETER  n=P**2
        if nargin<2:
            v_['LAMBDA'] = float(4.0);  #  SIF file default value
        else:
            v_['LAMBDA'] = float(args[1])
        v_['1.0'] = 1.0
        v_['1'] = 1
        v_['2'] = 2
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = v_['1.0']/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['C'] = v_['H2']*v_['LAMBDA']
        v_['-C'] = -1.0*v_['C']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2x_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['J+1'] = 1+J
                v_['J-1'] = -1+J
                [ig,ig_,_] = s2x_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
                iv = ix_['U'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
                iv = ix_['U'+str(int(v_['I+1']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['U'+str(int(v_['I-1']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['U'+str(I)+','+str(int(v_['J+1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['U'+str(I)+','+str(int(v_['J-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            pb.xlower[ix_['U'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['U'+str(int(v_['1']))+','+str(J)]] = 0.0
            pb.xlower[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
            pb.xupper[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['P']))]] = 0.0
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['P']))]] = 0.0
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2x_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
                    ielftype = arrset( ielftype,ie,iet_['eEXP'])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-C']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        pb.pbclass = "NOR2-MN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPU = np.exp(EV_[0])
        f_   = EXPU
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EXPU
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EXPU
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

