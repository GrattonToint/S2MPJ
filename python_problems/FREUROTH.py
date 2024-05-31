from s2xlib import *
class  FREUROTH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FREUROTH
#    *********
# 
#    The Freudentstein and Roth test problem
# 
#    Source: problem 2 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Toint#33, Buckley#24
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SUR2-AN-V-0"
# 
#    N is the number of variables
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FREUROTH'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'FREUROTH'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(4);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   5000           $-PARAMETER
        v_['NGS'] = -1+v_['N']
        v_['1'] = 1
        v_['2'] = 2
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
        for I in range(int(v_['1']),int(v_['NGS'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2x_ii('R'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-2.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2x_ii('S'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['X'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(-14.0)+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['NGS'])+1):
            v_['I+1'] = 1+I
            pbm.gconst = arrset(pbm.gconst,ig_['R'+str(I)],float(13.0))
            pbm.gconst = arrset(pbm.gconst,ig_['S'+str(I)],float(29.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.x0[ix_['X'+str(int(v_['1']))]] = float(0.5)
        pb.x0[ix_['X'+str(int(v_['2']))]] = float(-2.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eFRDRTH', iet_)
        elftv = loaset(elftv,it,0,'ELV')
        elftp = []
        elftp = loaset(elftp,it,0,'COEFF')
        elftp = loaset(elftp,it,1,'XCOEFF')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['NGS'])+1):
            v_['I+1'] = 1+I
            ename = 'A'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eFRDRTH')
                ielftype = arrset( ielftype,ie,iet_['eFRDRTH'])
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ELV')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='COEFF')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(5.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='XCOEFF')
            loaset(pbm.elpar,ie,posep[0],float(-1.0))
            ename = 'B'+str(I)
            [ie,ie_,newelt] = s2x_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eFRDRTH')
                ielftype = arrset( ielftype,ie,iet_['eFRDRTH'])
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='ELV')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='COEFF')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
            posep = find(elftp[ielftype[ie]],lambda x:x=='XCOEFF')
            loaset(pbm.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['NGS'])+1):
            ig = ig_['R'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['A'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['S'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eFRDRTH(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TWOC = pbm.elpar[iel_][0]+pbm.elpar[iel_][0]
        ELV2 = EV_[0]*EV_[0]
        XCELV = pbm.elpar[iel_][1]*EV_[0]
        f_   = (pbm.elpar[iel_][0]+XCELV)*ELV2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TWOC*EV_[0]+3.0*pbm.elpar[iel_][1]*ELV2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = TWOC+6.0*XCELV
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

