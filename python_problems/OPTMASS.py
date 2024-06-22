from s2mpjlib import *
class  OPTMASS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTMASS
#    *********
# 
#    A constrained optimal control problem
#    adapted from Gawande and Dunn
# 
#    The problem is that of a particle of unit mass moving on a
#    frictionless plane under the action of a controlling force whose
#    magnitude may not exceed unity. At time=0, the particle moves through
#    the origin of the plane in the direction of the positive x-axis with
#    speed SPEED.  The cost function incorporates two conflicting control
#    objectives, namely: maximization of the particle's final (at time=1)
#    distance from the origin and minimization of its final speed.  By
#    increasing the  value of the penalty constant PEN, more stress can be
#    placed on the latter objective.
# 
#    Gawande and Dunn originally use a starting point (in the control
#    only) that is much closer to the solution than the one chosen
#    here.
# 
#    Source:
#    M. Gawande and J. Dunn,
#    "A Projected Newton Method in a Cartesian Product of Balls",
#    JOTA 59(1): 59-69, 1988.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "QQR2-AN-V-V"
# 
#    Number of discretization steps in the time interval
#    The number of variables is 6 * (N + 2) -2 , 4 of which are fixed.
# 
#           Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n = 70    original value
# IE N                   100            $-PARAMETER n = 610
# IE N                   200            $-PARAMETER n = 1210
# IE N                   500            $-PARAMETER n = 3010
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'OPTMASS'

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
            v_['N'] = int(10);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   1000           $-PARAMETER n = 6010
# IE N                   5000           $-PARAMETER n = 30010
        v_['SPEED'] = 0.01
        v_['PEN'] = 0.335
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['N+1'] = 1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['1/N'] = 1.0/v_['RN']
        v_['-1/N'] = -1.0*v_['1/N']
        v_['1/N2'] = v_['1/N']*v_['1/N']
        v_['-1/2N2'] = -0.5*v_['1/N2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['2'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(I),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(J)+','+str(I))
                [iv,ix_,_] = s2mpj_ii('V'+str(J)+','+str(I),ix_)
                pb.xnames=arrset(pb.xnames,iv,'V'+str(J)+','+str(I))
                [iv,ix_,_] = s2mpj_ii('F'+str(J)+','+str(I),ix_)
                pb.xnames=arrset(pb.xnames,iv,'F'+str(J)+','+str(I))
        for J in range(int(v_['1']),int(v_['2'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(J)+','+str(int(v_['N+1'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(J)+','+str(int(v_['N+1'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(J)+','+str(int(v_['N+1'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(J)+','+str(int(v_['N+1'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('F',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N+1'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['2'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'A'+str(J)+','+str(I))
                iv = ix_['X'+str(J)+','+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(J)+','+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['V'+str(J)+','+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(v_['-1/N'])+pbm.A[ig,iv]
                iv = ix_['F'+str(J)+','+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(v_['-1/2N2'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('B'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'B'+str(J)+','+str(I))
                iv = ix_['V'+str(J)+','+str(I)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['V'+str(J)+','+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['F'+str(J)+','+str(int(v_['I-1']))]
                pbm.A[ig,iv] = float(v_['-1/N'])+pbm.A[ig,iv]
        for I in range(int(v_['0']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'C'+str(I))
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
        for I in range(int(v_['0']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['C'+str(I)],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['X'+str(int(v_['1']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['1']))+','+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['X'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['X'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        pb.xlower[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]] = v_['SPEED']
        pb.xupper[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]] = v_['SPEED']
        pb.xlower[ix_['V'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        pb.xupper[ix_['V'+str(int(v_['2']))+','+str(int(v_['0']))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.0))
        if('V'+str(int(v_['1']))+','+str(int(v_['0'])) in ix_):
            pb.x0[ix_['V'+str(int(v_['1']))+','+str(int(v_['0']))]] = float(v_['SPEED'])
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['V'+str(int(v_['1']))+','+str(int(v_['0']))]),float(v_['SPEED'])))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'O1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'X'+str(int(v_['1']))+','+str(int(v_['N+1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'O2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'X'+str(int(v_['2']))+','+str(int(v_['N+1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'O3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'V'+str(int(v_['1']))+','+str(int(v_['N+1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'O4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
        ielftype = arrset(ielftype, ie, iet_["eSQ"])
        vname = 'V'+str(int(v_['2']))+','+str(int(v_['N+1']))
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['2'])+1):
                ename = 'D'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'F'+str(J)+','+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,0.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['F']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['PEN']))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['PEN']))
        for I in range(int(v_['0']),int(v_['N'])+1):
            ig = ig_['C'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['1']))+','+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(int(v_['2']))+','+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           -0.04647
# LO SOLTN(100)          ???
# LO SOLTN(200)          ???
# LO SOLTN(500)          ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "QQR2-AN-V-V"
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

