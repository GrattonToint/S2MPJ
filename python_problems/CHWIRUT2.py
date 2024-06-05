from s2mpjlib import *
class  CHWIRUT2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHWIRUT2
#    *********
# 
#    NIST Data fitting problem CHWIRUT2 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = exp[-b1*x]/(b2+b3*x) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Chwirut, D., NIST (197?).  
#      Ultrasonic Reference Block Study. 
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "NOR2-MN-3-54"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHWIRUT2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'CHWIRUT2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 54
        v_['N'] = 3
        v_['1'] = 1
        v_['X1'] = 0.5
        v_['X2'] = 1.0
        v_['X3'] = 1.75
        v_['X4'] = 3.75
        v_['X5'] = 5.75
        v_['X6'] = 0.875
        v_['X7'] = 2.25
        v_['X8'] = 3.25
        v_['X9'] = 5.25
        v_['X10'] = 0.75
        v_['X11'] = 1.75
        v_['X12'] = 2.75
        v_['X13'] = 4.75
        v_['X14'] = 0.625
        v_['X15'] = 1.25
        v_['X16'] = 2.25
        v_['X17'] = 4.25
        v_['X18'] = 0.5
        v_['X19'] = 3.0
        v_['X20'] = 0.75
        v_['X21'] = 3.0
        v_['X22'] = 1.5
        v_['X23'] = 6.0
        v_['X24'] = 3.0
        v_['X25'] = 6.0
        v_['X26'] = 1.5
        v_['X27'] = 3.0
        v_['X28'] = 0.5
        v_['X29'] = 2.0
        v_['X30'] = 4.0
        v_['X31'] = 0.75
        v_['X32'] = 2.0
        v_['X33'] = 5.0
        v_['X34'] = 0.75
        v_['X35'] = 2.25
        v_['X36'] = 3.75
        v_['X37'] = 5.75
        v_['X38'] = 3.0
        v_['X39'] = 0.75
        v_['X40'] = 2.5
        v_['X41'] = 4.0
        v_['X42'] = 0.75
        v_['X43'] = 2.5
        v_['X44'] = 4.0
        v_['X45'] = 0.75
        v_['X46'] = 2.5
        v_['X47'] = 4.0
        v_['X48'] = 0.5
        v_['X49'] = 6.0
        v_['X50'] = 3.0
        v_['X51'] = 0.5
        v_['X52'] = 2.75
        v_['X53'] = 0.5
        v_['X54'] = 1.75
        v_['Y1'] = 92.9
        v_['Y2'] = 57.1
        v_['Y3'] = 31.05
        v_['Y4'] = 11.5875
        v_['Y5'] = 8.025
        v_['Y6'] = 63.6
        v_['Y7'] = 21.4
        v_['Y8'] = 14.25
        v_['Y9'] = 8.475
        v_['Y10'] = 63.8
        v_['Y11'] = 26.8
        v_['Y12'] = 16.4625
        v_['Y13'] = 7.125
        v_['Y14'] = 67.3
        v_['Y15'] = 41.0
        v_['Y16'] = 21.15
        v_['Y17'] = 8.175
        v_['Y18'] = 81.50
        v_['Y19'] = 13.12
        v_['Y20'] = 59.9
        v_['Y21'] = 14.62
        v_['Y22'] = 32.9
        v_['Y23'] = 5.44
        v_['Y24'] = 12.56
        v_['Y25'] = 5.44
        v_['Y26'] = 32.0
        v_['Y27'] = 13.95
        v_['Y28'] = 75.8
        v_['Y29'] = 20.0
        v_['Y30'] = 10.42
        v_['Y31'] = 59.5
        v_['Y32'] = 21.67
        v_['Y33'] = 8.55
        v_['Y34'] = 62.0
        v_['Y35'] = 20.2
        v_['Y36'] = 7.76
        v_['Y37'] = 3.75
        v_['Y38'] = 11.81
        v_['Y39'] = 54.7
        v_['Y40'] = 23.7
        v_['Y41'] = 11.55
        v_['Y42'] = 61.3
        v_['Y43'] = 17.7
        v_['Y44'] = 8.74
        v_['Y45'] = 59.2
        v_['Y46'] = 16.3
        v_['Y47'] = 8.62
        v_['Y48'] = 81.0
        v_['Y49'] = 4.87
        v_['Y50'] = 14.62
        v_['Y51'] = 81.7
        v_['Y52'] = 17.17
        v_['Y53'] = 81.3
        v_['Y54'] = 28.9
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'F'+str(I))
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
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('B1' in ix_):
            pb.x0[ix_['B1']] = float(0.1)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B1']),float(0.1)))
        if('B2' in ix_):
            pb.x0[ix_['B2']] = float(0.01)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B2']),float(0.01)))
        if('B3' in ix_):
            pb.x0[ix_['B3']] = float(0.02)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['B3']),float(0.02)))
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE16', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eE16')
            ielftype = arrset(ielftype, ie, iet_["eE16"])
            vname = 'B1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'B3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='X')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-MN-3-54"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE16(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E = np.exp(-EV_[0]*pbm.elpar[iel_][0])
        EX = E*pbm.elpar[iel_][0]
        EX2 = EX*pbm.elpar[iel_][0]
        V2PV3X = EV_[1]+EV_[2]*pbm.elpar[iel_][0]
        V2PV32 = V2PV3X*V2PV3X
        V2PV33 = V2PV3X*V2PV32
        f_   = E/V2PV3X
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -EX/V2PV3X
            g_[1] = -E/V2PV32
            g_[2] = -EX/V2PV32
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = EX2/V2PV3X
                H_[0,1] = EX/V2PV32
                H_[1,0] = H_[0,1]
                H_[0,2] = EX2/V2PV32
                H_[2,0] = H_[0,2]
                H_[1,1] = 2.0*E/V2PV33
                H_[1,2] = 2.0*EX/V2PV33
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EX2/V2PV33
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

