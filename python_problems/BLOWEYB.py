from s2mpjlib import *
class  BLOWEYB(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BLOWEYB
#    *********
# 
#    A nonconvex quadratic program proposed by 
#    James Blowey (University of Durham)
# 
#    Given function v(s) and f(s) = v(s) + A(inv) v(s), s in [0,1],
#    minimize 
# 
#         (u(s) - v(s))(trans) ( A + A(inv) ) u(s) - (u(s) - v(s))(trans)f(s)
# 
#    where 
# 
#       u(s) in [-1,1] and int[0,1] u(s) ds = int[0,1] v(s) ds
# 
#    and A is the 
# 
#       "- Laplacian with Neumann boundary conditions on a uniform mesh"
# 
#    The troublesome term A(inv) u(s) is replaced by the additional 
#    variable w(s) and the constraint A w(s) = u(s)
# 
#    The function v(s) is chosen to be 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BLOWEYB'

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
#    classification = "QLR2-MN-V-V"
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n = 22, m = 12
# IE N                   100            $-PARAMETER  n = 202, m = 102
# IE N                   1000           $-PARAMETER  n = 2002, m = 1002
# IE N                   2000           $-PARAMETER  n = 4002, m = 2002
        if nargin<1:
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   4000           $-PARAMETER  n = 8002, m = 4002
# IE N                   8000           $-PARAMETER  n = 16002, m = 8002
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['5'] = 5
        v_['9'] = 9
        v_['10'] = 10
        v_['ONE'] = 1.0
        v_['-ONE'] = -1.0
        v_['TWO'] = 2.0
        v_['-TWO'] = -2.0
        v_['RN'] = float(v_['N'])
        v_['N**2'] = v_['RN']*v_['RN']
        v_['N-1'] = -1+v_['N']
        v_['1/N**2'] = v_['ONE']/v_['N**2']
        v_['-1/N**2'] = v_['-ONE']/v_['N**2']
        v_['-2/N**2'] = 2.0*v_['-1/N**2']
        v_['2N**2'] = 2.0*v_['N**2']
        v_['-2N**2'] = -2.0*v_['N**2']
        v_['N/10'] = int(np.fix(v_['N']/v_['10']))
        v_['N/5'] = int(np.fix(v_['N']/v_['5']))
        v_['NA'] = v_['N/10']
        v_['A'] = float(v_['NA'])
        v_['A'] = v_['A']/v_['RN']
        v_['NA+1'] = 1+v_['NA']
        v_['NB'] = v_['N/5']*v_['2']
        v_['B'] = float(v_['NB'])
        v_['B'] = v_['B']/v_['RN']
        v_['NB+1'] = 1+v_['NB']
        v_['NC'] = v_['N/5']*v_['3']
        v_['C'] = float(v_['NC'])
        v_['C'] = v_['C']/v_['RN']
        v_['NC+1'] = 1+v_['NC']
        v_['ND'] = v_['N/10']*v_['9']
        v_['D'] = float(v_['ND'])
        v_['D'] = v_['D']/v_['RN']
        v_['ND+1'] = 1+v_['ND']
        v_['INT'] = v_['ONE']
        v_['INT'] = v_['INT']+v_['A']
        v_['INT'] = v_['INT']+v_['B']
        v_['INT'] = v_['INT']-v_['C']
        v_['INT'] = v_['INT']-v_['D']
        v_['INT'] = v_['INT']*v_['RN']
        for I in range(int(v_['0']),int(v_['NA'])+1):
            v_['V'+str(I)] = 1.0
        v_['STEP'] = v_['B']-v_['A']
        v_['STEP'] = v_['STEP']*v_['RN']
        v_['STEP'] = v_['TWO']/v_['STEP']
        for I in range(int(v_['NA+1']),int(v_['NB'])+1):
            v_['J'] = I-v_['NA']
            v_['RJ'] = float(v_['J'])
            v_['VAL'] = v_['RJ']*v_['STEP']
            v_['VAL'] = v_['ONE']-v_['VAL']
            v_['V'+str(I)] = v_['VAL']
        for I in range(int(v_['NB+1']),int(v_['NC'])+1):
            v_['V'+str(I)] = -1.0
        v_['STEP'] = v_['D']-v_['C']
        v_['STEP'] = v_['STEP']*v_['RN']
        v_['STEP'] = v_['TWO']/v_['STEP']
        for I in range(int(v_['NC+1']),int(v_['ND'])+1):
            v_['J'] = I-v_['NC']
            v_['RJ'] = float(v_['J'])
            v_['VAL'] = v_['RJ']*v_['STEP']
            v_['VAL'] = v_['-ONE']+v_['VAL']
            v_['V'+str(I)] = v_['VAL']
        for I in range(int(v_['ND']),int(v_['N'])+1):
            v_['V'+str(I)] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['VAL'] = v_['V'+str(I)]*v_['-1/N**2']
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['VAL'])+pbm.A[ig,iv]
            v_['VAL'] = v_['V'+str(I)]*v_['-2/N**2']
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['VAL'])+pbm.A[ig,iv]
        v_['VAL'] = v_['V'+str(int(v_['1']))]-v_['V'+str(int(v_['0']))]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(v_['VAL'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['VAL'] = -2.0*v_['V'+str(I)]
            v_['VAL'] = v_['VAL']+v_['V'+str(int(v_['I-1']))]
            v_['VAL'] = v_['VAL']+v_['V'+str(int(v_['I+1']))]
            [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(v_['VAL'])+pbm.A[ig,iv]
        v_['VAL'] = v_['V'+str(int(v_['N-1']))]-v_['V'+str(int(v_['N']))]
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['U'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['VAL'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('INT',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'INT')
        iv = ix_['U'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['0'])))
        iv = ix_['U'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['U'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['0'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['0'])))
        iv = ix_['W'+str(int(v_['0']))]
        pbm.A[ig,iv] = float(v_['-1/N**2'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            [ig,ig_,_] = s2mpj_ii('CON'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CON'+str(I))
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(2.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['U'+str(int(v_['I-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(v_['-1/N**2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('INT',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'INT')
            iv = ix_['U'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['N'])))
        iv = ix_['U'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        iv = ix_['U'+str(int(v_['N-1']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('CON'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON'+str(int(v_['N'])))
        iv = ix_['W'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(v_['-1/N**2'])+pbm.A[ig,iv]
        [ig,ig_,_] = s2mpj_ii('INT',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'INT')
        iv = ix_['U'+str(int(v_['N']))]
        pbm.A[ig,iv] = float(0.5)+pbm.A[ig,iv]
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
        pbm.gconst = arrset(pbm.gconst,ig_['INT'],float(v_['INT']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.xlower[ix_['U'+str(I)]] = -1.0
            pb.xupper[ix_['U'+str(I)]] = 1.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        for I in range(int(v_['0']),int(v_['N'])+1):
            pb.x0[ix_['U'+str(I)]] = float(v_['V'+str(I)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'Z')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            ename = 'C'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ename = 'D'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'O'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
            ielftype = arrset(ielftype, ie, iet_["ePROD"])
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'D'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        ename = 'D'+str(int(v_['N']))
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'U'+str(int(v_['N']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
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
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O'+str(int(v_['0']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-TWO']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['0']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['ONE']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-TWO']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['TWO']))
        ig = ig_['OBJ']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['N']))])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['ONE']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/N**2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -4.93517D+02   $ N = 10 
# XL SOLUTION            -5.30009D+03   $ N = 100
# XL SOLUTION            -5.33674D+04   $ N = 1000
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
        pb.pbclass = "QLR2-MN-V-V"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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

