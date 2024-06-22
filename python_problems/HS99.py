from s2mpjlib import *
class  HS99(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS99
#    *********
# 
#    Source: problem 99 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "OOR2-AN-7-2"
# 
#    Constants
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS99'

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
        v_['1'] = 1
        v_['2'] = 2
        v_['7'] = 7
        v_['8'] = 8
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['7'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        pbm.gscale = arrset(pbm.gscale,ig,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('Q8E',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'Q8E')
        [ig,ig_,_] = s2mpj_ii('S8E',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'S8E')
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
        pbm.gconst = arrset(pbm.gconst,ig_['Q8E'],float(100000.0))
        pbm.gconst = arrset(pbm.gconst,ig_['S8E'],float(1000.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),1.58)
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.5))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eR8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        [it,iet_,_] = s2mpj_ii( 'eQ8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        [it,iet_,_] = s2mpj_ii( 'eS8T', iet_)
        elftv = loaset(elftv,it,0,'X1')
        elftv = loaset(elftv,it,1,'X2')
        elftv = loaset(elftv,it,2,'X3')
        elftv = loaset(elftv,it,3,'X4')
        elftv = loaset(elftv,it,4,'X5')
        elftv = loaset(elftv,it,5,'X6')
        elftv = loaset(elftv,it,6,'X7')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        ename = 'R8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eR8T')
        ielftype = arrset(ielftype, ie, iet_["eR8T"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Q8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eQ8T')
        ielftype = arrset(ielftype, ie, iet_["eQ8T"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'S8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eS8T')
        ielftype = arrset(ielftype, ie, iet_["eS8T"])
        vname = 'X1'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X1')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X2'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X2')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X5'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X5')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X6'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X6')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X7'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,1.58,0.5)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X7')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        pbm.grftype = arrset(pbm.grftype,ig,'gL2')
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['R8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['Q8E']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Q8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['S8E']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -831079892.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-7-2"
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
        pbm.efpar = arrset( pbm.efpar,0,50.0)
        pbm.efpar = arrset( pbm.efpar,1,50.0)
        pbm.efpar = arrset( pbm.efpar,2,75.0)
        pbm.efpar = arrset( pbm.efpar,3,75.0)
        pbm.efpar = arrset( pbm.efpar,4,75.0)
        pbm.efpar = arrset( pbm.efpar,5,100.0)
        pbm.efpar = arrset( pbm.efpar,6,100.0)
        pbm.efpar = arrset( pbm.efpar,7,25.0)
        pbm.efpar = arrset( pbm.efpar,8,25.0)
        pbm.efpar = arrset( pbm.efpar,9,50.0)
        pbm.efpar = arrset( pbm.efpar,10,50.0)
        pbm.efpar = arrset( pbm.efpar,11,50.0)
        pbm.efpar = arrset( pbm.efpar,12,90.0)
        pbm.efpar = arrset( pbm.efpar,13,90.0)
        pbm.efpar = arrset( pbm.efpar,14,32.0)
        return pbm

    @staticmethod
    def eR8T(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R2 = pbm.efpar[0]*pbm.efpar[7]*np.cos(EV_[0])
        R3 = pbm.efpar[1]*pbm.efpar[8]*np.cos(EV_[1])+R2
        R4 = pbm.efpar[2]*pbm.efpar[9]*np.cos(EV_[2])+R3
        R5 = pbm.efpar[3]*pbm.efpar[10]*np.cos(EV_[3])+R4
        R6 = pbm.efpar[4]*pbm.efpar[11]*np.cos(EV_[4])+R5
        R7 = pbm.efpar[5]*pbm.efpar[12]*np.cos(EV_[5])+R6
        f_   = pbm.efpar[6]*pbm.efpar[13]*np.cos(EV_[6])+R7
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -pbm.efpar[0]*pbm.efpar[7]*np.sin(EV_[0])
            g_[1] = -pbm.efpar[1]*pbm.efpar[8]*np.sin(EV_[1])
            g_[2] = -pbm.efpar[2]*pbm.efpar[9]*np.sin(EV_[2])
            g_[3] = -pbm.efpar[3]*pbm.efpar[10]*np.sin(EV_[3])
            g_[4] = -pbm.efpar[4]*pbm.efpar[11]*np.sin(EV_[4])
            g_[5] = -pbm.efpar[5]*pbm.efpar[12]*np.sin(EV_[5])
            g_[6] = -pbm.efpar[6]*pbm.efpar[13]*np.sin(EV_[6])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = -pbm.efpar[0]*pbm.efpar[7]*np.cos(EV_[0])
                H_[1,1] = -pbm.efpar[1]*pbm.efpar[8]*np.cos(EV_[1])
                H_[2,2] = -pbm.efpar[2]*pbm.efpar[9]*np.cos(EV_[2])
                H_[3,3] = -pbm.efpar[3]*pbm.efpar[10]*np.cos(EV_[3])
                H_[4,4] = -pbm.efpar[4]*pbm.efpar[11]*np.cos(EV_[4])
                H_[5,5] = -pbm.efpar[5]*pbm.efpar[12]*np.cos(EV_[5])
                H_[6,6] = -pbm.efpar[6]*pbm.efpar[13]*np.cos(EV_[6])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eS8T(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S2 = pbm.efpar[7]*(pbm.efpar[0]*np.sin(EV_[0])-pbm.efpar[14])
        S3 = pbm.efpar[8]*(pbm.efpar[1]*np.sin(EV_[1])-pbm.efpar[14])+S2
        S4 = pbm.efpar[9]*(pbm.efpar[2]*np.sin(EV_[2])-pbm.efpar[14])+S3
        S5 = pbm.efpar[10]*(pbm.efpar[3]*np.sin(EV_[3])-pbm.efpar[14])+S4
        S6 = pbm.efpar[11]*(pbm.efpar[4]*np.sin(EV_[4])-pbm.efpar[14])+S5
        S7 = pbm.efpar[12]*(pbm.efpar[5]*np.sin(EV_[5])-pbm.efpar[14])+S6
        f_   = pbm.efpar[13]*(pbm.efpar[6]*np.sin(EV_[6])-pbm.efpar[14])+S7
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.efpar[0]*pbm.efpar[7]*np.cos(EV_[0])
            g_[1] = pbm.efpar[1]*pbm.efpar[8]*np.cos(EV_[1])
            g_[2] = pbm.efpar[2]*pbm.efpar[9]*np.cos(EV_[2])
            g_[3] = pbm.efpar[3]*pbm.efpar[10]*np.cos(EV_[3])
            g_[4] = pbm.efpar[4]*pbm.efpar[11]*np.cos(EV_[4])
            g_[5] = pbm.efpar[5]*pbm.efpar[12]*np.cos(EV_[5])
            g_[6] = pbm.efpar[6]*pbm.efpar[13]*np.cos(EV_[6])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = -pbm.efpar[0]*pbm.efpar[7]*np.sin(EV_[0])
                H_[1,1] = -pbm.efpar[1]*pbm.efpar[8]*np.sin(EV_[1])
                H_[2,2] = -pbm.efpar[2]*pbm.efpar[9]*np.sin(EV_[2])
                H_[3,3] = -pbm.efpar[3]*pbm.efpar[10]*np.sin(EV_[3])
                H_[4,4] = -pbm.efpar[4]*pbm.efpar[11]*np.sin(EV_[4])
                H_[5,5] = -pbm.efpar[5]*pbm.efpar[12]*np.sin(EV_[5])
                H_[6,6] = -pbm.efpar[6]*pbm.efpar[13]*np.sin(EV_[6])
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQ8T(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S2 = pbm.efpar[7]*(pbm.efpar[0]*np.sin(EV_[0])-pbm.efpar[14])
        S3 = pbm.efpar[8]*(pbm.efpar[1]*np.sin(EV_[1])-pbm.efpar[14])+S2
        S4 = pbm.efpar[9]*(pbm.efpar[2]*np.sin(EV_[2])-pbm.efpar[14])+S3
        S5 = pbm.efpar[10]*(pbm.efpar[3]*np.sin(EV_[3])-pbm.efpar[14])+S4
        S6 = pbm.efpar[11]*(pbm.efpar[4]*np.sin(EV_[4])-pbm.efpar[14])+S5
        S7 = pbm.efpar[12]*(pbm.efpar[5]*np.sin(EV_[5])-pbm.efpar[14])+S6
        DSD1 = pbm.efpar[0]*pbm.efpar[7]*np.cos(EV_[0])
        DSD2 = pbm.efpar[1]*pbm.efpar[8]*np.cos(EV_[1])
        DSD3 = pbm.efpar[2]*pbm.efpar[9]*np.cos(EV_[2])
        DSD4 = pbm.efpar[3]*pbm.efpar[10]*np.cos(EV_[3])
        DSD5 = pbm.efpar[4]*pbm.efpar[11]*np.cos(EV_[4])
        DSD6 = pbm.efpar[5]*pbm.efpar[12]*np.cos(EV_[5])
        DSD7 = pbm.efpar[6]*pbm.efpar[13]*np.cos(EV_[6])
        D2SD1 = -pbm.efpar[0]*pbm.efpar[7]*np.sin(EV_[0])
        D2SD2 = -pbm.efpar[1]*pbm.efpar[8]*np.sin(EV_[1])
        D2SD3 = -pbm.efpar[2]*pbm.efpar[9]*np.sin(EV_[2])
        D2SD4 = -pbm.efpar[3]*pbm.efpar[10]*np.sin(EV_[3])
        D2SD5 = -pbm.efpar[4]*pbm.efpar[11]*np.sin(EV_[4])
        D2SD6 = -pbm.efpar[5]*pbm.efpar[12]*np.sin(EV_[5])
        D2SD7 = -pbm.efpar[6]*pbm.efpar[13]*np.sin(EV_[6])
        Q2  = (
              0.5*pbm.efpar[7]*pbm.efpar[7]*(pbm.efpar[0]*np.sin(EV_[0])-pbm.efpar[14]))
        Q3 = (0.5*pbm.efpar[8]*pbm.efpar[8]*(pbm.efpar[1]*np.sin(EV_[1])-pbm.efpar[14])+
             pbm.efpar[8]*S2+Q2)
        Q4 = (0.5*pbm.efpar[9]*pbm.efpar[9]*(pbm.efpar[2]*np.sin(EV_[2])-pbm.efpar[14])+
             pbm.efpar[9]*S3+Q3)
        Q5  = (
              0.5*pbm.efpar[10]*pbm.efpar[10]*(pbm.efpar[3]*np.sin(EV_[3])-pbm.efpar[14])+pbm.efpar[10]*S4+Q4)
        Q6  = (
              0.5*pbm.efpar[11]*pbm.efpar[11]*(pbm.efpar[4]*np.sin(EV_[4])-pbm.efpar[14])+pbm.efpar[11]*S5+Q5)
        Q7  = (
              0.5*pbm.efpar[12]*pbm.efpar[12]*(pbm.efpar[5]*np.sin(EV_[5])-pbm.efpar[14])+pbm.efpar[12]*S6+Q6)
        f_    = (
              0.5*pbm.efpar[13]*pbm.efpar[13]*(pbm.efpar[6]*np.sin(EV_[6])-pbm.efpar[14])+pbm.efpar[13]*S7+Q7)
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (0.5*pbm.efpar[7]*pbm.efpar[7]*pbm.efpar[0]*np.cos(EV_[0])+
                 (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9]+pbm.efpar[8])*DSD1)
            g_[1] = (0.5*pbm.efpar[8]*pbm.efpar[8]*pbm.efpar[1]*np.cos(EV_[1])+
                 (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9])*DSD2)
            g_[2] = (0.5*pbm.efpar[9]*pbm.efpar[9]*pbm.efpar[2]*np.cos(EV_[2])+
                 (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10])*DSD3)
            g_[3] = (0.5*pbm.efpar[10]*pbm.efpar[10]*pbm.efpar[3]*np.cos(EV_[3])+
                 (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11])*DSD4)
            g_[4] = (0.5*pbm.efpar[11]*pbm.efpar[11]*pbm.efpar[4]*np.cos(EV_[4])+
                 (pbm.efpar[13]+pbm.efpar[12])*DSD5)
            g_[5] = (0.5*pbm.efpar[12]*pbm.efpar[12]*pbm.efpar[5]*np.cos(EV_[5])+
                 pbm.efpar[13]*DSD6)
            g_[6] = 0.5*pbm.efpar[13]*pbm.efpar[13]*pbm.efpar[6]*np.cos(EV_[6])
            if nargout>2:
                H_ = np.zeros((7,7))
                H_[0,0] = (-0.5*pbm.efpar[7]*pbm.efpar[7]*pbm.efpar[0]*np.sin(EV_[0])+
                     (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9]+pbm.efpar[8])*D2SD1)
                H_[1,1] = (-0.5*pbm.efpar[8]*pbm.efpar[8]*pbm.efpar[1]*np.sin(EV_[1])+
                     (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10]+pbm.efpar[9])*D2SD2)
                H_[2,2] = (-0.5*pbm.efpar[9]*pbm.efpar[9]*pbm.efpar[2]*np.sin(EV_[2])+
                     (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11]+pbm.efpar[10])*D2SD3)
                H_[3,3] = (-0.5*pbm.efpar[10]*pbm.efpar[10]*pbm.efpar[3]*np.sin(EV_[3])+
                     (pbm.efpar[13]+pbm.efpar[12]+pbm.efpar[11])*D2SD4)
                H_[4,4] = (-0.5*pbm.efpar[11]*pbm.efpar[11]*pbm.efpar[4]*np.sin(EV_[4])+
                     (pbm.efpar[13]+pbm.efpar[12])*D2SD5)
                H_[5,5] = (-0.5*pbm.efpar[12]*pbm.efpar[12]*pbm.efpar[5]*np.sin(EV_[5])+
                     pbm.efpar[13]*D2SD6)
                H_[6,6] = -0.5*pbm.efpar[13]*pbm.efpar[13]*pbm.efpar[6]*np.sin(EV_[6])
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

