from s2xlib import *
class  PENLT2NE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PENLT2NE
#    --------
# 
#    The second penalty function
# 
#    This is a nonlinear least-squares problem with M=2*N groups.
#     Group 1 is linear.
#     Groups 2 to N use 2 nonlinear elements.
#     Groups N+1 to M-1 use 1 nonlinear element.
#     Group M uses N nonlinear elements.
#    The Hessian matrix is dense. This is a nonlinear equation version
#    of problem PENALTY2.
# 
#    Source:  Problem 24 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#112 (p. 80)
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "NOR2-AN-V-V"
# 
#    Number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PENLT2NE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'PENLT2NE'
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
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
        v_['A'] = 0.00001
        v_['B'] = 1.0
        v_['1'] = 1
        v_['2'] = 2
        v_['N+1'] = 1+v_['N']
        v_['M'] = v_['N']+v_['N']
        v_['M-1'] = -1+v_['M']
        v_['EM1/10'] = np.exp(-0.1)
        v_['1/A'] = 1.0/v_['A']
        v_['1/B'] = 1.0/v_['B']
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
        [ig,ig_,_] = s2x_ii('G'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G'+str(int(v_['1'])))
        iv = ix_['X'+str(int(v_['1']))]
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        [ig,ig_,_] = s2x_ii('G'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G'+str(int(v_['1'])))
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/B']))
        for I in range(int(v_['2']),int(v_['M-1'])+1):
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(I))
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/A']))
        [ig,ig_,_] = s2x_ii('G'+str(int(v_['M'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'G'+str(int(v_['M'])))
        pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/B']))
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
        pbm.gconst = arrset(pbm.gconst,ig_['G1'],float(0.2))
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['RI-1'] = float(v_['I-1'])
            v_['I/10'] = 0.1*v_['RI']
            v_['I-1/10'] = 0.1*v_['RI-1']
            v_['EI/10'] = np.exp(v_['I/10'])
            v_['EI-1/10'] = np.exp(v_['I-1/10'])
            v_['YI'] = v_['EI/10']+v_['EI-1/10']
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['YI']))
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['EM1/10']))
        pbm.gconst = arrset(pbm.gconst,ig_['G'+str(int(v_['M']))],float(1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE10', iet_)
        elftv = loaset(elftv,it,0,'V')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['2']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ename = 'A'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE10')
            ielftype = arrset(ielftype, ie, iet_["eE10"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE10')
            ielftype = arrset(ielftype, ie, iet_["eE10"])
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            v_['-N'] = -1*v_['N']
            v_['I-N'] = I+v_['-N']
            v_['I-N+1'] = 1+v_['I-N']
            ename = 'C'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE10')
            ielftype = arrset(ielftype, ie, iet_["eE10"])
            vname = 'X'+str(int(v_['I-N+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['N'])+1):
            ename = 'D'+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(J)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['N'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        for I in range(int(v_['N+1']),int(v_['M-1'])+1):
            ig = ig_['G'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for J in range(int(v_['1']),int(v_['N'])+1):
            v_['-J'] = -1*J
            v_['N-J'] = v_['N']+v_['-J']
            v_['N-J+1'] = 1+v_['N-J']
            v_['WI'] = float(v_['N-J+1'])
            ig = ig_['G'+str(int(v_['M']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['D'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['WI']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
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
        pb.pbclass = "NOR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE10(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EXPA = np.exp(0.1*EV_[0])
        f_   = EXPA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.1*EXPA
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.01*EXPA
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

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

