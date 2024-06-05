from s2mpjlib import *
class  HS111(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS111
#    *********
# 
#    This problem is a chemical equilibrium problem involving 3 linear
#    equality constraints.
# 
#    Source: problem 111 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "OOR2-AN-10-3"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS111'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS111'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 10
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
        v_['C1'] = -6.089
        v_['C2'] = -17.164
        v_['C3'] = -34.054
        v_['C4'] = -5.914
        v_['C5'] = -24.721
        v_['C6'] = -14.986
        v_['C7'] = -24.100
        v_['C8'] = -10.708
        v_['C9'] = -26.662
        v_['C10'] = -22.179
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
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
        [ig,ig_,_] = s2mpj_ii('CON1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON1')
        [ig,ig_,_] = s2mpj_ii('CON2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON2')
        [ig,ig_,_] = s2mpj_ii('CON3',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CON3')
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
        pbm.gconst = arrset(pbm.gconst,ig_['CON1'],float(2.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON2'],float(1.0))
        pbm.gconst = arrset(pbm.gconst,ig_['CON3'],float(1.0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-100.0)
        pb.xupper = np.full((pb.n,1),100.0)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(-2.3))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOBJ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        elftp = []
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'eEXP', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'+str(I)] = float(I)
            v_['RI'+str(I)] = 0.1+v_['RI'+str(I)]
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'O'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eOBJ')
            ielftype = arrset(ielftype, ie, iet_["eOBJ"])
            v_['TEMP'] = v_['RI'+str(int(v_['1']))]
            v_['RI'+str(int(v_['1']))] = v_['RI'+str(I)]
            v_['RI'+str(I)] = v_['TEMP']
            v_['R'] = v_['RI'+str(int(v_['1']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['2']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['3']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['4']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V4')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['5']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V5')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['6']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V6')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['7']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V7')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['8']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V8')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['9']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V9')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            v_['R'] = v_['RI'+str(int(v_['10']))]
            v_['J'] = int(np.fix(v_['R']))
            vname = 'X'+str(int(v_['J']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V10')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            posep = find(elftp[ielftype[ie]],lambda x:x=='C')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(I)]))
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEXP')
            ielftype = arrset(ielftype, ie, iet_["eEXP"])
            vname = 'X'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,-100.0,100.0,-2.3)
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
            ig = ig_['OBJ']
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['O'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        ig = ig_['CON1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E2'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        ig = ig_['CON2']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E4'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E5'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E6'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        ig = ig_['CON3']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E3'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E7'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E8'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        posel = posel+1
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E9'])
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(2.0))
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E10'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -47.707579
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-AN-10-3"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        EX = np.exp(EV_[0])
        f_   = EX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EX
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = EX
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOBJ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E1 = np.exp(EV_[0])
        E2 = np.exp(EV_[1])
        E3 = np.exp(EV_[2])
        E4 = np.exp(EV_[3])
        E5 = np.exp(EV_[4])
        E6 = np.exp(EV_[5])
        E7 = np.exp(EV_[6])
        E8 = np.exp(EV_[7])
        E9 = np.exp(EV_[8])
        E10 = np.exp(EV_[9])
        SUM = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10
        f_   = E1*(pbm.elpar[iel_][0]+EV_[0]-np.log(SUM))
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E1*(pbm.elpar[iel_][0]+EV_[0]-np.log(SUM))+E1*(1.0e+0-E1/SUM)
            g_[1] = -E1*E2/SUM
            g_[2] = -E1*E3/SUM
            g_[3] = -E1*E4/SUM
            g_[4] = -E1*E5/SUM
            g_[5] = -E1*E6/SUM
            g_[6] = -E1*E7/SUM
            g_[7] = -E1*E8/SUM
            g_[8] = -E1*E9/SUM
            g_[9] = -E1*E10/SUM
            if nargout>2:
                H_ = np.zeros((10,10))
                H_[0,0] = (E1*(pbm.elpar[iel_][0]+EV_[0]-np.log(SUM))+E1*(1.0e+0-E1/SUM)+
                     E1*(1.0e+0-E1/SUM)+E1*(-E1/SUM)+E1*(E1**2/SUM**2))
                H_[0,1] = (-1.0e+0+E1/SUM)*E1*E2/SUM
                H_[1,0] = H_[0,1]
                H_[1,1] = (-1.0e+0+E2/SUM)*E1*E2/SUM
                H_[0,2] = (-1.0e+0+E1/SUM)*E1*E3/SUM
                H_[2,0] = H_[0,2]
                H_[1,2] = E1*E2*E3/SUM**2
                H_[2,1] = H_[1,2]
                H_[2,2] = (-1.0e+0+E3/SUM)*E1*E3/SUM
                H_[0,3] = (-1.0e+0+E1/SUM)*E1*E4/SUM
                H_[3,0] = H_[0,3]
                H_[1,3] = E1*E2*E4/SUM**2
                H_[3,1] = H_[1,3]
                H_[2,3] = E1*E3*E4/SUM**2
                H_[3,2] = H_[2,3]
                H_[3,3] = (-1.0e+0+E4/SUM)*E1*E4/SUM
                H_[0,4] = (-1.0e+0+E1/SUM)*E1*E5/SUM
                H_[4,0] = H_[0,4]
                H_[1,4] = E1*E2*E5/SUM**2
                H_[4,1] = H_[1,4]
                H_[2,4] = E1*E3*E5/SUM**2
                H_[4,2] = H_[2,4]
                H_[3,4] = E1*E4*E5/SUM**2
                H_[4,3] = H_[3,4]
                H_[4,4] = (-1.0e+0+E5/SUM)*E1*E5/SUM
                H_[0,5] = (-1.0e+0+E1/SUM)*E1*E6/SUM
                H_[5,0] = H_[0,5]
                H_[1,5] = E1*E2*E6/SUM**2
                H_[5,1] = H_[1,5]
                H_[2,5] = E1*E3*E6/SUM**2
                H_[5,2] = H_[2,5]
                H_[3,5] = E1*E4*E6/SUM**2
                H_[5,3] = H_[3,5]
                H_[4,5] = E1*E5*E6/SUM**2
                H_[5,4] = H_[4,5]
                H_[5,5] = (-1.0e+0+E6/SUM)*E1*E6/SUM
                H_[0,6] = (-1.0e+0+E1/SUM)*E1*E7/SUM
                H_[6,0] = H_[0,6]
                H_[1,6] = E1*E2*E7/SUM**2
                H_[6,1] = H_[1,6]
                H_[2,6] = E1*E3*E7/SUM**2
                H_[6,2] = H_[2,6]
                H_[3,6] = E1*E4*E7/SUM**2
                H_[6,3] = H_[3,6]
                H_[4,6] = E1*E5*E7/SUM**2
                H_[6,4] = H_[4,6]
                H_[5,6] = E1*E6*E7/SUM**2
                H_[6,5] = H_[5,6]
                H_[6,6] = (-1.0e+0+E7/SUM)*E1*E7/SUM
                H_[0,7] = (-1.0e+0+E1/SUM)*E1*E8/SUM
                H_[7,0] = H_[0,7]
                H_[1,7] = E1*E2*E8/SUM**2
                H_[7,1] = H_[1,7]
                H_[2,7] = E1*E3*E8/SUM**2
                H_[7,2] = H_[2,7]
                H_[3,7] = E1*E4*E8/SUM**2
                H_[7,3] = H_[3,7]
                H_[4,7] = E1*E5*E8/SUM**2
                H_[7,4] = H_[4,7]
                H_[5,7] = E1*E6*E8/SUM**2
                H_[7,5] = H_[5,7]
                H_[6,7] = E1*E7*E8/SUM**2
                H_[7,6] = H_[6,7]
                H_[7,7] = (-1.0e+0+E8/SUM)*E1*E8/SUM
                H_[0,8] = (-1.0e+0+E1/SUM)*E1*E9/SUM
                H_[8,0] = H_[0,8]
                H_[1,8] = E1*E2*E9/SUM**2
                H_[8,1] = H_[1,8]
                H_[2,8] = E1*E3*E9/SUM**2
                H_[8,2] = H_[2,8]
                H_[3,8] = E1*E4*E9/SUM**2
                H_[8,3] = H_[3,8]
                H_[4,8] = E1*E5*E9/SUM**2
                H_[8,4] = H_[4,8]
                H_[5,8] = E1*E6*E9/SUM**2
                H_[8,5] = H_[5,8]
                H_[6,8] = E1*E7*E9/SUM**2
                H_[8,6] = H_[6,8]
                H_[7,8] = E1*E8*E9/SUM**2
                H_[8,7] = H_[7,8]
                H_[8,8] = (-1.0e+0+E9/SUM)*E1*E9/SUM
                H_[0,9] = (-1.0e+0+E1/SUM)*E1*E10/SUM
                H_[9,0] = H_[0,9]
                H_[1,9] = E1*E2*E10/SUM**2
                H_[9,1] = H_[1,9]
                H_[2,9] = E1*E3*E10/SUM**2
                H_[9,2] = H_[2,9]
                H_[3,9] = E1*E4*E10/SUM**2
                H_[9,3] = H_[3,9]
                H_[4,9] = E1*E5*E10/SUM**2
                H_[9,4] = H_[4,9]
                H_[5,9] = E1*E6*E10/SUM**2
                H_[9,5] = H_[5,9]
                H_[6,9] = E1*E7*E10/SUM**2
                H_[9,6] = H_[6,9]
                H_[7,9] = E1*E8*E10/SUM**2
                H_[9,7] = H_[7,9]
                H_[8,9] = E1*E9*E10/SUM**2
                H_[9,8] = H_[8,9]
                H_[9,9] = (-1.0e+0+E10/SUM)*E1*E10/SUM
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

