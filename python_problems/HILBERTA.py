from s2xlib import *
class  HILBERTA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HILBERTA
#    *********
# 
#    The Hilbert quadratic
# 
#    Source:
#    K. Schittkowski,
#    "More Test Examples for Nonlinear Programming Codes",
#    Springer Verlag, Heidelberg, 1987.
# 
#    See also Buckley#19 (p. 59)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "QUR2-AN-V-0"
# 
#    Dimension of the problem
# 
#           Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER Schittkowski 274
# IE N                   4              $-PARAMETER Schittkowski 275
# IE N                   5              $-PARAMETER Buckley 19
# IE N                   6              $-PARAMETER Schittkowski 276
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HILBERTA'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HILBERTA'
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
        if nargin<2:
            v_['D'] = float(0.0);  #  SIF file default value
        else:
            v_['D'] = float(args[1])
        v_['1'] = 1
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(I)+1):
                [ig,ig_,_] = s2x_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(-3.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PR')
                ielftype = arrset(ielftype, ie, iet_["en2PR"])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-3.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-3.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'E'+str(I)+','+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
            ielftype = arrset(ielftype, ie, iet_["eSQ"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,-3.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['I+J'] = I+J
                v_['I+J-1'] = -1+v_['I+J']
                v_['RINVH'] = float(v_['I+J-1'])
                v_['HIJ'] = 1.0/v_['RINVH']
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['HIJ']))
            v_['2I'] = 2*I
            v_['2I-1'] = -1+v_['2I']
            v_['RH0'] = float(v_['2I-1'])
            v_['HII'] = 1.0/v_['RH0']
            v_['HII/2'] = 0.5*v_['HII']
            v_['COEFF'] = v_['HII/2']+v_['D']
            ig = ig_['G'+str(I)+','+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['COEFF']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "QUR2-AN-V-0"
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

