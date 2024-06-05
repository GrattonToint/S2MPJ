from s2mpjlib import *
class  CHEMRCTA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHEMRCTA
#    *********
# 
#    The tubular chemical reactor model problem by Poore, using a
#    finite difference approximation to the steady state solutions.
# 
#    Source: Problem 8, eqs (8.6)--(8.9) in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               minor correction by Ph. Shott, Jan 1995 and F Ruediger, Mar 1997.
# 
#    classification = "NOR2-MN-V-V"
# 
#    The axial coordinate interval is [0,1]
# 
#    Number of discretized point for the interval [0,1].
#    The number of variables is 2N.
# 
#           Alternative values for the SIF file parameters:
# IE N                   5              $-PARAMETER n = 10
# IE N                   25             $-PARAMETER n = 50
# IE N                   50             $-PARAMETER n = 100
# IE N                   250            $-PARAMETER n = 500    original value
# IE N                   500            $-PARAMETER n = 1000
# IE N                   2500           $-PARAMETER n = 5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHEMRCTA'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CHEMRCTA'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        if nargin<2:
            v_['PEM'] = float(1.0);  #  SIF file default value
        else:
            v_['PEM'] = float(args[1])
        if nargin<3:
            v_['PEH'] = float(5.0);  #  SIF file default value
        else:
            v_['PEH'] = float(args[2])
        if nargin<4:
            v_['D'] = float(0.135);  #  SIF file default value
        else:
            v_['D'] = float(args[3])
        if nargin<5:
            v_['B'] = float(0.5);  #  SIF file default value
        else:
            v_['B'] = float(args[4])
        if nargin<6:
            v_['BETA'] = float(2.0);  #  SIF file default value
        else:
            v_['BETA'] = float(args[5])
        if nargin<7:
            v_['GAMMA'] = float(25.0);  #  SIF file default value
        else:
            v_['GAMMA'] = float(args[6])
        v_['1'] = 1
        v_['2'] = 2
        v_['1.0'] = 1.0
        v_['N-1'] = -1+v_['N']
        v_['1/H'] = float(v_['N-1'])
        v_['-1/H'] = -1.0*v_['1/H']
        v_['H'] = v_['1.0']/v_['1/H']
        v_['1/H2'] = v_['1/H']*v_['1/H']
        v_['-D'] = -1.0*v_['D']
        v_['1/PEM'] = v_['1.0']/v_['PEM']
        v_['1/H2PEM'] = v_['1/PEM']*v_['1/H2']
        v_['-1/H2PM'] = -1.0*v_['1/H2PEM']
        v_['HPEM'] = v_['PEM']*v_['H']
        v_['-HPEM'] = -1.0*v_['HPEM']
        v_['-2/H2PM'] = v_['-1/H2PM']+v_['-1/H2PM']
        v_['CU1'] = 1.0*v_['-HPEM']
        v_['CUI-1'] = v_['1/H2PEM']+v_['1/H']
        v_['CUI'] = v_['-2/H2PM']+v_['-1/H']
        v_['BD'] = v_['B']*v_['D']
        v_['-BETA'] = -1.0*v_['BETA']
        v_['1/PEH'] = v_['1.0']/v_['PEH']
        v_['1/H2PEH'] = v_['1/PEH']*v_['1/H2']
        v_['-1/H2PH'] = -1.0*v_['1/H2PEH']
        v_['HPEH'] = v_['PEH']*v_['H']
        v_['-HPEH'] = -1.0*v_['HPEH']
        v_['-2/H2PH'] = v_['-1/H2PH']+v_['-1/H2PH']
        v_['CT1'] = 1.0*v_['-HPEH']
        v_['CTI-1'] = v_['1/H2PEH']+v_['1/H']
        v_['CTI'] = v_['-2/H2PH']+v_['-1/H']
        v_['CTI'] = v_['CTI']+v_['-BETA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'T'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['1'])))
        iv = ix_['U'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['1'])))
        iv = ix_['U'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(v_['CU1'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['1'])))
        iv = ix_['T'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['1'])))
        iv = ix_['T'+str(int(v_['2']))]
        pbm.A[ig,iv] = float(v_['CT1'])+pbm.A[ig,iv]
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('GU'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'GU'+str(I))
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['CUI-1'])+pbm.A[ig,iv]
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['CUI'])+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['1/H2PEM'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('GT'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'GT'+str(I))
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(v_['BETA'])+pbm.A[ig,iv]
            iv = ix_['T'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(v_['CTI-1'])+pbm.A[ig,iv]
            iv = ix_['T'+str(I)]
            pbm.A[ig,iv] = float(v_['CTI'])+pbm.A[ig,iv]
            iv = ix_['T'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(v_['1/H2PEH'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['N'])))
        iv = ix_['U'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['N'])))
        iv = ix_['U'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['N'])))
        iv = ix_['T'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['N'])))
        iv = ix_['T'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['GU'+str(int(v_['1']))],float(v_['-HPEM'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['GT'+str(int(v_['1']))],float(v_['-HPEH'])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['T'+str(I)]] = 0.0000001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eREAC', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'T')
        elftp = []
        elftp = loaset(elftp,it,0,'G')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            ename = 'EU'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eREAC')
                ielftype = arrset( ielftype,ie,iet_['eREAC'])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['GAMMA']))
            ename = 'ET'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eREAC')
                ielftype = arrset( ielftype,ie,iet_['eREAC'])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'T'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='G')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['GAMMA']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            ig = ig_['GU'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EU'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-D']))
            ig = ig_['GT'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['ET'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['BD']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eREAC(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DADT = pbm.elpar[iel_][0]/(EV_[1]*EV_[1])
        D2ADT2 = -2.0*DADT/EV_[1]
        EX = np.exp(pbm.elpar[iel_][0]-pbm.elpar[iel_][0]/EV_[1])
        UEX = EX*EV_[0]
        f_   = UEX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EX
            g_[1] = UEX*DADT
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EX*DADT
                H_[1,0] = H_[0,1]
                H_[1,1] = UEX*(DADT*DADT+D2ADT2)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

