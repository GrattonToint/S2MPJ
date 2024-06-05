from s2mpjlib import *
class  POROUS2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POROUS2
#    *********
# 
#    The problem is to solve the porous medium equation on the unit square.
#    The equation is
# 
#        \Delta ( u^2 ) + d \frac{\partial}{\partial x_1}( u^3 ) + f = 0
# 
#    within the domain.  The boundary condition are that u = 1 on the bottom
#    and left sides and u = 0 on the top and right sides.  Discretization is
#    using the usual central differences. The function f is a point source of
#    maginitude 50 at the lower left grid point.  The initial approximation
#    is a discretization of 1 - x_1 x_2.
# 
#    Source: example 3.2.4 in
#    S. Eisenstat and H. Walker,
#    "Choosing the forcing terms in an inexact Newton method"
#    Report 6/94/75, Dept of Maths, Utah State University, 1994.
# 
#    SIF input: Ph. Toint, July 1994.
# 
#    classification = "NOR2-MN-V-V"
# 
#    P is the number of points in one side of the unit square.
#    There are P*P variables.
# 
#           Alternative values for the SIF file parameters:
# IE P                   32             $-PARAMETER      original value
# IE P                   64             $-PARAMETER 
# IE P                   72             $-PARAMETER 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'POROUS2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'POROUS2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['P'] = int(32);  #  SIF file default value
        else:
            v_['P'] = int(args[0])
        if nargin<2:
            v_['D'] = float(-50.0);  #  SIF file default value
        else:
            v_['D'] = float(args[1])
        v_['1'] = 1
        v_['2'] = 2
        v_['P-1'] = -1+v_['P']
        v_['RP-1'] = float(v_['P-1'])
        v_['H'] = 1.0/v_['RP-1']
        v_['H2'] = v_['H']*v_['H']
        v_['1/H2'] = 1.0/v_['H2']
        v_['2H'] = 2.0*v_['H']
        v_['D/2H'] = v_['D']/v_['2H']
        v_['-D/2H'] = -1.0*v_['D/2H']
        v_['-4/H2'] = -4.0*v_['1/H2']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['P'])+1):
            for I in range(int(v_['1']),int(v_['P'])+1):
                [iv,ix_,_] = s2mpj_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
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
              arrset(pbm.gconst,ig_['G'+str(int(v_['P-1']))+','+str(int(v_['P-1']))],float(-50.0)))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for J in range(int(v_['1']),int(v_['P'])+1):
            pb.xlower[ix_['U'+str(int(v_['1']))+','+str(J)]] = 1.0
            pb.xupper[ix_['U'+str(int(v_['1']))+','+str(J)]] = 1.0
            pb.xlower[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
            pb.xupper[ix_['U'+str(int(v_['P']))+','+str(J)]] = 0.0
        for I in range(int(v_['2']),int(v_['P-1'])+1):
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['P']))]] = 1.0
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['P']))]] = 1.0
            pb.xlower[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
            pb.xupper[ix_['U'+str(I)+','+str(int(v_['1']))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['RI'] = float(I)
                v_['RJ'] = float(J)
                v_['I-1'] = -1.0+v_['RI']
                v_['J-1'] = -1.0+v_['RJ']
                v_['X1'] = v_['I-1']*v_['H']
                v_['X2'] = v_['J-1']*v_['H']
                v_['X1X2'] = v_['X1']*v_['X2']
                v_['MX1X2'] = -1.0*v_['X1X2']
                v_['UIJ'] = 1.0+v_['MX1X2']
                pb.x0[ix_['U'+str(I)+','+str(J)]] = float(v_['UIJ'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        [it,iet_,_] = s2mpj_ii( 'eCB', iet_)
        elftv = loaset(elftv,it,0,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                ename = 'US'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eSQ')
                ielftype = arrset(ielftype, ie, iet_["eSQ"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'UC'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eCB')
                ielftype = arrset(ielftype, ie, iet_["eCB"])
                vname = 'U'+str(I)+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
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
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            for J in range(int(v_['2']),int(v_['P-1'])+1):
                v_['J-1'] = -1+J
                v_['J+1'] = 1+J
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['US'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['US'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['US'+str(I)+','+str(int(v_['J-1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['US'+str(I)+','+str(int(v_['J+1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['US'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-4/H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['UC'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['D/2H']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['UC'+str(int(v_['I-1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-D/2H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
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
    def eCB(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 3.0*EV_[0]*EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 6.0*EV_[0]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

