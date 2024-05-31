from s2xlib import *
class  SPARSINE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPARSINE
#    *********
# 
#    A sparse problem involving sine functions
# 
#    SIF input: Nick Gould, November 1995
# 
#    classification = "OUR2-AN-V-0"
# 
#    The number of variables 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPARSINE'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'SPARSINE'
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
#           Alternative values for the SIF file parameters:
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
        v_['1'] = 1
        v_['2'] = 2
        v_['4'] = 4
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
            [ig,ig_,_] = s2x_ii('OBJ'+str(I),ig_)
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
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eSINE', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eSINE')
            ielftype = arrset(ielftype, ie, iet_["eSINE"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.5)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gSQR',igt_)
        [it,igt_,_] = s2x_ii('gSQR',igt_)
        grftp = []
        grftp = loaset(grftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        pbm.grpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['OBJ'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gSQR')
            v_['RI'] = float(I)
            posgp = find(grftp[igt_[pbm.grftype[ig]]],lambda x:x=='P')
            pbm.grpar =loaset(pbm.grpar,ig,posgp[0],float(v_['RI']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            v_['J'] = 2*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            v_['J'] = 3*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            v_['J'] = 5*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            v_['J'] = 7*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
            v_['J'] = 11*I
            v_['J'] = -1+v_['J']
            v_['K'] = int(np.fix(v_['J']/v_['N']))
            v_['K'] = v_['K']*v_['N']
            v_['J'] = v_['J']-v_['K']
            v_['J'] = 1+v_['J']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(int(v_['J']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSINE(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = np.sin(EV_[0])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(EV_[0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -np.sin(EV_[0])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQR(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 0.5*pbm.grpar[igr_][0]*GVAR_*GVAR_
        if nargout>1:
            g_ = pbm.grpar[igr_][0]*GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = pbm.grpar[igr_][0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

