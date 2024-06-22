from s2mpjlib import *
class  RAYBENDL(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A ray bending problem.  A ray across a inhomogeneous 2D medium is
#    represented by a piecewise linear curve whose knots can be chosen.  
#    The problem is then to optimize the positions of these knots in order 
#    to obtain a ray path corresponding to the minimum travel time from 
#    source to receiver,  according to Fermat principle.
# 
#    The problem becomes harder and harder when the dimesnion increases
#    because the knots are getting closer and closer and the objective
#    has a nondifferentiable kink when two knots coincide.  The difficulty
#    is less apparent when exact second derivatives are not used.
# 
#    Source: a test example in
#    T.J. Moser, G. Nolet and R. Snieder,
#    "Ray bending revisited",
#    Bulletin of the Seism. Society of America 21(1).
# 
#    SIF input: Ph Toint, Dec 1991.
# 
#    classification = "OXR2-MY-V-0"
# 
#    number of  knots  ( >= 4 )
#    ( n = 2( NKNOTS - 1 ) ) 
# 
#           Alternative values for the SIF file parameters:
# IE NKNOTS              4              $-PARAMETER n = 6
# IE NKNOTS              11             $-PARAMETER n = 20
# IE NKNOTS              21             $-PARAMETER n = 40     original value
# IE NKNOTS              32             $-PARAMETER n = 62
# IE NKNOTS              64             $-PARAMETER n = 126
# IE NKNOTS              512            $-PARAMETER n = 1022
# IE NKNOTS              1024           $-PARAMETER n = 2046
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'RAYBENDL'

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
            v_['NKNOTS'] = int(4);  #  SIF file default value
        else:
            v_['NKNOTS'] = int(args[0])
        v_['XSRC'] = 0.0
        v_['ZSRC'] = 0.0
        v_['XRCV'] = 100.0
        v_['ZRCV'] = 100.0
        v_['NK-1'] = -1+v_['NKNOTS']
        v_['NK-2'] = -2+v_['NKNOTS']
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['NKNOTS'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2mpj_ii('Z'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Z'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['NKNOTS'])+1):
            [ig,ig_,_] = s2mpj_ii('TIME'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            pbm.gscale = arrset(pbm.gscale,ig,float(2.0))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        pb.xlower[ix_['X'+str(int(v_['0']))]] = v_['XSRC']
        pb.xupper[ix_['X'+str(int(v_['0']))]] = v_['XSRC']
        pb.xlower[ix_['Z'+str(int(v_['0']))]] = v_['ZSRC']
        pb.xupper[ix_['Z'+str(int(v_['0']))]] = v_['ZSRC']
        pb.xlower[ix_['X'+str(int(v_['NKNOTS']))]] = v_['XRCV']
        pb.xupper[ix_['X'+str(int(v_['NKNOTS']))]] = v_['XRCV']
        pb.xlower[ix_['Z'+str(int(v_['NKNOTS']))]] = v_['ZRCV']
        pb.xupper[ix_['Z'+str(int(v_['NKNOTS']))]] = v_['ZRCV']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        v_['XRANGE'] = v_['XRCV']-v_['XSRC']
        v_['ZRANGE'] = v_['ZRCV']-v_['ZSRC']
        v_['RKNOTS'] = float(v_['NKNOTS'])
        for I in range(int(v_['0']),int(v_['NKNOTS'])+1):
            v_['REALI'] = float(I)
            v_['FRAC'] = v_['REALI']/v_['RKNOTS']
            v_['XINCR'] = v_['FRAC']*v_['XRANGE']
            v_['ZINCR'] = v_['FRAC']*v_['ZRANGE']
            v_['XC'] = v_['XSRC']+v_['XINCR']
            v_['ZC'] = v_['ZSRC']+v_['ZINCR']
            pb.x0[ix_['X'+str(I)]] = float(v_['XC'])
            pb.x0[ix_['Z'+str(I)]] = float(v_['ZC'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eTT', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'Z1')
        elftv = loaset(elftv,it,3,'Z2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NKNOTS'])+1):
            v_['I-1'] = -1+I
            ename = 'T'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                pbm.elftype = arrset(pbm.elftype,ie,'eTT')
                ielftype = arrset( ielftype,ie,iet_['eTT'])
            vname = 'X'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z'+str(int(v_['I-1']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'Z'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Z2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NKNOTS'])+1):
            ig = ig_['TIME'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['T'+str(I)])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#   Solution of the continuous problem
# LO RAYBENDL            96.2424
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OXR2-MY-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,0.01)
        return pbm

    @staticmethod
    def eTT(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,4))
        IV_ = np.zeros(3)
        U_[2,0] = U_[2,0]-1
        U_[2,1] = U_[2,1]+1
        U_[0,2] = U_[0,2]+1
        U_[1,3] = U_[1,3]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        IV_[2] = U_[2:3,:].dot(EV_)
        C0 = 1.0+pbm.efpar[0]*IV_[0]
        C1 = 1.0+pbm.efpar[0]*IV_[1]
        DCDZ = pbm.efpar[0]
        V = 1.0/C1+1.0/C0
        VDZ0 = -DCDZ/(C0*C0)
        VDZ1 = -DCDZ/(C1*C1)
        VDZ0Z0 = 2.0*DCDZ*DCDZ/C0**3
        VDZ1Z1 = 2.0*DCDZ*DCDZ/C1**3
        DZ1 = IV_[1]-IV_[0]
        R = np.sqrt(IV_[2]*IV_[2]+DZ1*DZ1)
        RDX = IV_[2]/R
        RDZ1 = DZ1/R
        RDZ0 = -RDZ1
        RDXDX = (1.0-IV_[2]*IV_[2]/(R*R))/R
        RDXZ1 = -IV_[2]*DZ1/R**3
        RDXZ0 = -RDXZ1
        RDZ1Z1 = (1.0-DZ1*DZ1/(R*R))/R
        RDZ0Z0 = RDZ1Z1
        RDZ0Z1 = -RDZ1Z1
        f_   = V*R
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[2] = V*RDX
            g_[0] = V*RDZ0+VDZ0*R
            g_[1] = V*RDZ1+VDZ1*R
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[2,2] = V*RDXDX
                H_[2,0] = VDZ0*RDX+V*RDXZ0
                H_[0,2] = H_[2,0]
                H_[2,1] = VDZ1*RDX+V*RDXZ1
                H_[1,2] = H_[2,1]
                H_[0,0] = V*RDZ0Z0+VDZ0Z0*R+2.0*VDZ0*RDZ0
                H_[0,1] = V*RDZ0Z1+VDZ1*RDZ0+VDZ0*RDZ1
                H_[1,0] = H_[0,1]
                H_[1,1] = V*RDZ1Z1+VDZ1Z1*R+2.0*VDZ1*RDZ1
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

