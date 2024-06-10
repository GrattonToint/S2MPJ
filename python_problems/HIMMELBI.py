from s2mpjlib import *
class  HIMMELBI(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBI
#    *********
# 
#    An unpleasant weapon assignment problem by Bracken and McCormick.
# 
#    The real problem has integer variables.
#    Also, the sign of ci have been reversed in order to have a
#    meaningful constraints on the total number of weapons (a fully
#    desirable situation).
# 
#    Source: problem 23 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "OLR2-MN-100-12"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HIMMELBI'

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
        v_['N'] = 100
        v_['NT'] = 20
        v_['A1,1'] = 1.0
        v_['A2,1'] = 0.84
        v_['A3,1'] = 0.96
        v_['A4,1'] = 1.0
        v_['A5,1'] = 0.92
        v_['A1,2'] = 0.95
        v_['A2,2'] = 0.83
        v_['A3,2'] = 0.95
        v_['A4,2'] = 1.0
        v_['A5,2'] = 0.94
        v_['A1,3'] = 1.0
        v_['A2,3'] = 0.85
        v_['A3,3'] = 0.95
        v_['A4,3'] = 1.0
        v_['A5,3'] = 0.92
        v_['A1,4'] = 1.0
        v_['A2,4'] = 0.84
        v_['A3,4'] = 0.96
        v_['A4,4'] = 1.0
        v_['A5,4'] = 0.95
        v_['A1,5'] = 1.0
        v_['A2,5'] = 0.85
        v_['A3,5'] = 0.96
        v_['A4,5'] = 1.0
        v_['A5,5'] = 0.95
        v_['A1,6'] = 0.85
        v_['A2,6'] = 0.81
        v_['A3,6'] = 0.90
        v_['A4,6'] = 1.0
        v_['A5,6'] = 0.98
        v_['A1,7'] = 0.90
        v_['A2,7'] = 0.81
        v_['A3,7'] = 0.92
        v_['A4,7'] = 1.0
        v_['A5,7'] = 0.98
        v_['A1,8'] = 0.85
        v_['A2,8'] = 0.82
        v_['A3,8'] = 0.91
        v_['A4,8'] = 1.0
        v_['A5,8'] = 1.0
        v_['A1,9'] = 0.80
        v_['A2,9'] = 0.80
        v_['A3,9'] = 0.92
        v_['A4,9'] = 1.0
        v_['A5,9'] = 1.0
        v_['A1,10'] = 1.0
        v_['A2,10'] = 0.86
        v_['A3,10'] = 0.95
        v_['A4,10'] = 0.96
        v_['A5,10'] = 0.90
        v_['A1,11'] = 1.0
        v_['A2,11'] = 1.0
        v_['A3,11'] = 0.99
        v_['A4,11'] = 0.91
        v_['A5,11'] = 0.95
        v_['A1,12'] = 1.0
        v_['A2,12'] = 0.98
        v_['A3,12'] = 0.98
        v_['A4,12'] = 0.92
        v_['A5,12'] = 0.96
        v_['A1,13'] = 1.0
        v_['A2,13'] = 1.0
        v_['A3,13'] = 0.99
        v_['A4,13'] = 0.91
        v_['A5,13'] = 0.91
        v_['A1,14'] = 1.0
        v_['A2,14'] = 0.88
        v_['A3,14'] = 0.98
        v_['A4,14'] = 0.92
        v_['A5,14'] = 0.98
        v_['A1,15'] = 1.0
        v_['A2,15'] = 0.87
        v_['A3,15'] = 0.97
        v_['A4,15'] = 0.98
        v_['A5,15'] = 0.99
        v_['A1,16'] = 1.0
        v_['A2,16'] = 0.88
        v_['A3,16'] = 0.98
        v_['A4,16'] = 0.93
        v_['A5,16'] = 0.99
        v_['A1,17'] = 1.0
        v_['A2,17'] = 0.85
        v_['A3,17'] = 0.95
        v_['A4,17'] = 1.0
        v_['A5,17'] = 1.0
        v_['A1,18'] = 0.95
        v_['A2,18'] = 0.84
        v_['A3,18'] = 0.92
        v_['A4,18'] = 1.0
        v_['A5,18'] = 1.0
        v_['A1,19'] = 1.0
        v_['A2,19'] = 0.85
        v_['A3,19'] = 0.93
        v_['A4,19'] = 1.0
        v_['A5,19'] = 1.0
        v_['A1,20'] = 1.0
        v_['A2,20'] = 0.85
        v_['A3,20'] = 0.92
        v_['A4,20'] = 1.0
        v_['A5,20'] = 1.0
        v_['B1'] = 30.0
        v_['B6'] = 100.0
        v_['B10'] = 40.0
        v_['B14'] = 50.0
        v_['B15'] = 70.0
        v_['B16'] = 35.0
        v_['B20'] = 10.0
        v_['U1'] = 60.0
        v_['U2'] = 50.0
        v_['U3'] = 50.0
        v_['U4'] = 75.0
        v_['U5'] = 40.0
        v_['U6'] = 60.0
        v_['U7'] = 35.0
        v_['U8'] = 30.0
        v_['U9'] = 25.0
        v_['U10'] = 150.0
        v_['U11'] = 30.0
        v_['U12'] = 45.0
        v_['U13'] = 125.0
        v_['U14'] = 200.0
        v_['U15'] = 200.0
        v_['U16'] = 130.0
        v_['U17'] = 100.0
        v_['U18'] = 100.0
        v_['U19'] = 100.0
        v_['U20'] = 150.0
        v_['C1'] = 200.0
        v_['C2'] = 100.0
        v_['C3'] = 300.0
        v_['C4'] = 150.0
        v_['C5'] = 250.0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['NW'] = 0.0
        for I in range(int(v_['1']),int(v_['5'])+1):
            v_['NW'] = v_['NW']+v_['C'+str(I)]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['1']),int(v_['NT'])+1):
            for I in range(int(v_['1']),int(v_['5'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['1']),int(v_['NT'])+1):
            [ig,ig_,_] = s2mpj_ii('P'+str(J),ig_)
            gtype = arrset(gtype,ig,'<>')
            v_['1/UJ'] = 1.0/v_['U'+str(J)]
            pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/UJ']))
        v_['J'] = 1
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 6
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 10
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 14
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 15
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 16
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        v_['J'] = 20
        for I in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('CB'+str(int(v_['J'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'CB'+str(int(v_['J'])))
            iv = ix_['X'+str(I)+','+str(int(v_['J']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['5'])+1):
            for J in range(int(v_['1']),int(v_['NT'])+1):
                [ig,ig_,_] = s2mpj_ii('CC'+str(I),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'CC'+str(I))
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        for J in range(int(v_['1']),int(v_['NT'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['P'+str(J)],float(1.0))
        v_['J'] = 1
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 6
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 10
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 14
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 15
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 16
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        v_['J'] = 20
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['CB'+str(int(v_['J']))],float(v_['B'+str(int(v_['J']))])))
        for I in range(int(v_['1']),int(v_['5'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CC'+str(I)],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),v_['NW'])
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en5PEXP', iet_)
        elftv = loaset(elftv,it,0,'Y1')
        elftv = loaset(elftv,it,1,'Y2')
        elftv = loaset(elftv,it,2,'Y3')
        elftv = loaset(elftv,it,3,'Y4')
        elftv = loaset(elftv,it,4,'Y5')
        elftp = []
        elftp = loaset(elftp,it,0,'A1')
        elftp = loaset(elftp,it,1,'A2')
        elftp = loaset(elftp,it,2,'A3')
        elftp = loaset(elftp,it,3,'A4')
        elftp = loaset(elftp,it,4,'A5')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for J in range(int(v_['1']),int(v_['NT'])+1):
            ename = 'PP'+str(J)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en5PEXP')
            ielftype = arrset(ielftype, ie, iet_["en5PEXP"])
            pb.x0 = np.zeros((pb.n,1))
            vname = 'X'+str(int(v_['1']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,v_['NW'],None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['2']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,v_['NW'],None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['3']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,v_['NW'],None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['4']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,v_['NW'],None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['5']))+','+str(J)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,v_['NW'],None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='Y5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='A1')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(int(v_['1']))+','+str(J)])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='A2')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(int(v_['2']))+','+str(J)])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='A3')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(int(v_['3']))+','+str(J)])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='A4')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(int(v_['4']))+','+str(J)])))
            posep = find(elftp[ielftype[ie]],lambda x:x=='A5')
            pbm.elpar  = (
                  loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(int(v_['5']))+','+str(J)])))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['1']),int(v_['NT'])+1):
            ig = ig_['P'+str(J)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PP'+str(J)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1735.56958
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OLR2-MN-100-12"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en5PEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        LA1 = np.log(pbm.elpar[iel_][0])
        LA2 = np.log(pbm.elpar[iel_][1])
        LA3 = np.log(pbm.elpar[iel_][2])
        LA4 = np.log(pbm.elpar[iel_][3])
        LA5 = np.log(pbm.elpar[iel_][4])
        FF  = (
              pbm.elpar[iel_][0]**EV_[0]*pbm.elpar[iel_][1]**EV_[1]*pbm.elpar[iel_][2]**EV_[2]*pbm.elpar[iel_][3]**EV_[3]*pbm.elpar[iel_][4]**EV_[4])
        f_   = FF
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = LA1*FF
            g_[1] = LA2*FF
            g_[2] = LA3*FF
            g_[3] = LA4*FF
            g_[4] = LA5*FF
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = LA1*LA1*FF
                H_[0,1] = LA1*LA2*FF
                H_[1,0] = H_[0,1]
                H_[0,2] = LA1*LA3*FF
                H_[2,0] = H_[0,2]
                H_[0,3] = LA1*LA4*FF
                H_[3,0] = H_[0,3]
                H_[0,4] = LA1*LA5*FF
                H_[4,0] = H_[0,4]
                H_[1,1] = LA2*LA2*FF
                H_[1,2] = LA2*LA3*FF
                H_[2,1] = H_[1,2]
                H_[1,3] = LA2*LA4*FF
                H_[3,1] = H_[1,3]
                H_[1,4] = LA2*LA5*FF
                H_[4,1] = H_[1,4]
                H_[2,2] = LA3*LA3*FF
                H_[2,3] = LA3*LA4*FF
                H_[3,2] = H_[2,3]
                H_[2,4] = LA3*LA5*FF
                H_[4,2] = H_[2,4]
                H_[3,3] = LA4*LA4*FF
                H_[3,4] = LA4*LA5*FF
                H_[4,3] = H_[3,4]
                H_[4,4] = LA5*LA5*FF
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

