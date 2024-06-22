from s2mpjlib import *
class  HS119(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS119
#    *********
# 
#    Source: problem 119 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Original Source: problem 7 in
#    A.R. Colville
#    "A comparative study on nonlinear programming"
#    IBM Scientific Center Report 320-2949, New York, 1968.
# 
#    SIF input: A.R. Conn, March 1991.
# 
#    classification = "OLR2-AN-16-8"
# 
#    Set useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS119'

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
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['15'] = 15
        v_['16'] = 16
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                v_['A'+str(I)+','+str(J)] = 0.0
        for I in range(int(v_['1']),int(v_['8'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                v_['B'+str(I)+','+str(J)] = 0.0
        for I in range(int(v_['1']),int(v_['16'])+1):
            v_['A'+str(I)+','+str(I)] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['4']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['8']))] = 1.0
        v_['A'+str(int(v_['1']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['3']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['2']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['9']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['3']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['7']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['11']))] = 1.0
        v_['A'+str(int(v_['4']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['6']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['12']))] = 1.0
        v_['A'+str(int(v_['5']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['6']))+','+str(int(v_['8']))] = 1.0
        v_['A'+str(int(v_['6']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['7']))+','+str(int(v_['11']))] = 1.0
        v_['A'+str(int(v_['7']))+','+str(int(v_['13']))] = 1.0
        v_['A'+str(int(v_['8']))+','+str(int(v_['10']))] = 1.0
        v_['A'+str(int(v_['8']))+','+str(int(v_['15']))] = 1.0
        v_['A'+str(int(v_['9']))+','+str(int(v_['12']))] = 1.0
        v_['A'+str(int(v_['9']))+','+str(int(v_['16']))] = 1.0
        v_['A'+str(int(v_['10']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['11']))+','+str(int(v_['13']))] = 1.0
        v_['A'+str(int(v_['12']))+','+str(int(v_['14']))] = 1.0
        v_['A'+str(int(v_['13']))+','+str(int(v_['14']))] = 1.0
        v_['B'+str(int(v_['1']))+','+str(int(v_['1']))] = 0.22
        v_['B'+str(int(v_['2']))+','+str(int(v_['1']))] = -1.46
        v_['B'+str(int(v_['3']))+','+str(int(v_['1']))] = 1.29
        v_['B'+str(int(v_['4']))+','+str(int(v_['1']))] = -1.10
        v_['B'+str(int(v_['7']))+','+str(int(v_['1']))] = 1.12
        v_['B'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.20
        v_['B'+str(int(v_['3']))+','+str(int(v_['2']))] = -0.89
        v_['B'+str(int(v_['4']))+','+str(int(v_['2']))] = -1.06
        v_['B'+str(int(v_['6']))+','+str(int(v_['2']))] = -1.72
        v_['B'+str(int(v_['8']))+','+str(int(v_['2']))] = 0.45
        v_['B'+str(int(v_['1']))+','+str(int(v_['3']))] = 0.19
        v_['B'+str(int(v_['2']))+','+str(int(v_['3']))] = -1.30
        v_['B'+str(int(v_['4']))+','+str(int(v_['3']))] = 0.95
        v_['B'+str(int(v_['6']))+','+str(int(v_['3']))] = -0.33
        v_['B'+str(int(v_['8']))+','+str(int(v_['3']))] = 0.26
        v_['B'+str(int(v_['1']))+','+str(int(v_['4']))] = 0.25
        v_['B'+str(int(v_['2']))+','+str(int(v_['4']))] = 1.82
        v_['B'+str(int(v_['4']))+','+str(int(v_['4']))] = -0.54
        v_['B'+str(int(v_['5']))+','+str(int(v_['4']))] = -1.43
        v_['B'+str(int(v_['7']))+','+str(int(v_['4']))] = 0.31
        v_['B'+str(int(v_['8']))+','+str(int(v_['4']))] = -1.10
        v_['B'+str(int(v_['1']))+','+str(int(v_['5']))] = 0.15
        v_['B'+str(int(v_['2']))+','+str(int(v_['5']))] = -1.15
        v_['B'+str(int(v_['3']))+','+str(int(v_['5']))] = -1.16
        v_['B'+str(int(v_['5']))+','+str(int(v_['5']))] = 1.51
        v_['B'+str(int(v_['6']))+','+str(int(v_['5']))] = 1.62
        v_['B'+str(int(v_['8']))+','+str(int(v_['5']))] = 0.58
        v_['B'+str(int(v_['1']))+','+str(int(v_['6']))] = 0.11
        v_['B'+str(int(v_['3']))+','+str(int(v_['6']))] = -0.96
        v_['B'+str(int(v_['4']))+','+str(int(v_['6']))] = -1.78
        v_['B'+str(int(v_['5']))+','+str(int(v_['6']))] = 0.59
        v_['B'+str(int(v_['6']))+','+str(int(v_['6']))] = 1.24
        v_['B'+str(int(v_['1']))+','+str(int(v_['7']))] = 0.12
        v_['B'+str(int(v_['2']))+','+str(int(v_['7']))] = 0.80
        v_['B'+str(int(v_['4']))+','+str(int(v_['7']))] = -0.41
        v_['B'+str(int(v_['5']))+','+str(int(v_['7']))] = -0.33
        v_['B'+str(int(v_['6']))+','+str(int(v_['7']))] = 0.21
        v_['B'+str(int(v_['7']))+','+str(int(v_['7']))] = 1.12
        v_['B'+str(int(v_['8']))+','+str(int(v_['7']))] = -1.03
        v_['B'+str(int(v_['1']))+','+str(int(v_['8']))] = 0.13
        v_['B'+str(int(v_['3']))+','+str(int(v_['8']))] = -0.49
        v_['B'+str(int(v_['5']))+','+str(int(v_['8']))] = -0.43
        v_['B'+str(int(v_['6']))+','+str(int(v_['8']))] = -0.26
        v_['B'+str(int(v_['8']))+','+str(int(v_['8']))] = 0.10
        v_['B'+str(int(v_['1']))+','+str(int(v_['9']))] = 1.00
        v_['B'+str(int(v_['7']))+','+str(int(v_['9']))] = -0.36
        v_['B'+str(int(v_['2']))+','+str(int(v_['10']))] = 1.00
        v_['B'+str(int(v_['3']))+','+str(int(v_['11']))] = 1.00
        v_['B'+str(int(v_['4']))+','+str(int(v_['12']))] = 1.00
        v_['B'+str(int(v_['5']))+','+str(int(v_['13']))] = 1.00
        v_['B'+str(int(v_['6']))+','+str(int(v_['14']))] = 1.00
        v_['B'+str(int(v_['7']))+','+str(int(v_['15']))] = 1.00
        v_['B'+str(int(v_['8']))+','+str(int(v_['16']))] = 1.00
        v_['C'+str(int(v_['1']))] = 2.5
        v_['C'+str(int(v_['2']))] = 1.1
        v_['C'+str(int(v_['3']))] = -3.1
        v_['C'+str(int(v_['4']))] = -3.5
        v_['C'+str(int(v_['5']))] = 1.3
        v_['C'+str(int(v_['6']))] = 2.1
        v_['C'+str(int(v_['7']))] = 2.3
        v_['C'+str(int(v_['8']))] = -1.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            [ig,ig_,_] = s2mpj_ii('OG'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['8'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(I))
                iv = ix_['X'+str(J)]
                pbm.A[ig,iv] = float(v_['B'+str(I)+','+str(J)])+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['8'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['C'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = np.full((pb.n,1),5.0)
        pb.xlower = np.zeros((pb.n,1))
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(10.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'U1')
        elftv = loaset(elftv,it,1,'U2')
        elftp = []
        elftp = loaset(elftp,it,0,'AIJ')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                ename = 'S'+str(I)+','+str(J)
                [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
                if newelt:
                    pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                    ielftype = arrset( ielftype,ie,iet_['ePROD'])
                vname = 'X'+str(I)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,5.0,10.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,5.0,10.0)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='AIJ')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['A'+str(I)+','+str(J)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['16'])+1):
            for J in range(int(v_['1']),int(v_['16'])+1):
                ig = ig_['OG'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['S'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
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
        pb.pbclass = "OLR2-AN-16-8"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TU1P1 = 2.0*EV_[0]+1
        TU2P1 = 2.0*EV_[1]+1
        FIRST = EV_[0]**2+EV_[0]+1.0
        SECOND = EV_[1]**2+EV_[1]+1.0
        f_   = pbm.elpar[iel_][0]*FIRST*SECOND
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*TU1P1*SECOND
            g_[1] = pbm.elpar[iel_][0]*TU2P1*FIRST
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = pbm.elpar[iel_][0]*2.0*SECOND
                H_[0,1] = pbm.elpar[iel_][0]*TU1P1*TU2P1
                H_[1,0] = H_[0,1]
                H_[1,1] = pbm.elpar[iel_][0]*2.0*FIRST
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

