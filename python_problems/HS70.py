from s2mpjlib import *
class  HS70(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS70
#    *********
# 
#    This problem arises in water flow routing.
# 
#    Source: problem 70 incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991, modified May 2024
# 
#    classification = "SQR2-MN-4-1"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'HS70'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'HS70'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['N'] = 4
        v_['1'] = 1
        v_['2'] = 2
        v_['19'] = 19
        v_['C1'] = 0.1
        v_['C2'] = 1.0
        v_['C3'] = 2.0
        v_['C4'] = 3.0
        v_['C5'] = 4.0
        v_['C6'] = 5.0
        v_['C7'] = 6.0
        v_['C8'] = 7.0
        v_['C9'] = 8.0
        v_['C10'] = 9.0
        v_['C11'] = 10.0
        v_['C12'] = 11.0
        v_['C13'] = 12.0
        v_['C14'] = 13.0
        v_['C15'] = 14.0
        v_['C16'] = 15.0
        v_['C17'] = 16.0
        v_['C18'] = 17.0
        v_['C19'] = 18.0
        v_['Y1'] = 0.00189
        v_['Y2'] = 0.1038
        v_['Y3'] = 0.268
        v_['Y4'] = 0.506
        v_['Y5'] = 0.577
        v_['Y6'] = 0.604
        v_['Y7'] = 0.725
        v_['Y8'] = 0.898
        v_['Y9'] = 0.947
        v_['Y10'] = 0.845
        v_['Y11'] = 0.702
        v_['Y12'] = 0.528
        v_['Y13'] = 0.385
        v_['Y14'] = 0.257
        v_['Y15'] = 0.159
        v_['Y16'] = 0.0869
        v_['Y17'] = 0.0453
        v_['Y18'] = 0.01509
        v_['Y19'] = 0.00189
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
        for I in range(int(v_['1']),int(v_['19'])+1):
            [ig,ig_,_] = s2mpj_ii('OBJ'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('C1',ig_)
        gtype = arrset(gtype,ig,'>=')
        cnames = arrset(cnames,ig,'C1')
        iv = ix_['X3']
        pbm.A[ig,iv] = float(1.0e+0)+pbm.A[ig,iv]
        iv = ix_['X4']
        pbm.A[ig,iv] = float(1.0e+0)+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['19'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['OBJ'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.00001)
        pb.xupper = np.full((pb.n,1),100.0)
        pb.xupper[ix_['X3']] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(2.0)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(2.0)))
        if('X2' in ix_):
            pb.x0[ix_['X2']] = float(4.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X2']),float(4.0))
        if('X3' in ix_):
            pb.x0[ix_['X3']] = float(0.04)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X3']),float(0.04)))
        if('X4' in ix_):
            pb.x0[ix_['X4']] = float(2.0)
        else:
            pb.y0 = arrset(pb.y0,find(pbm.congrps,lambda x:x==ig_['X4']),float(2.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eY1', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = []
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'eY2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'C')
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'X3')
        elftv = loaset(elftv,it,1,'X4')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['19'])+1):
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eY1')
            ielftype = arrset(ielftype, ie, iet_["eY1"])
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X2'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='C')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(I)]))
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eY2')
            ielftype = arrset(ielftype, ie, iet_["eY2"])
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X1'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X3'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'X4'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Y'+str(I)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='C')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['C'+str(I)]))
        ename = 'C1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
        ielftype = arrset(ielftype, ie, iet_["ePROD"])
        vname = 'X3'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X3')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'X4'
        [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,0.00001,100.0,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X4')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQR',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['19'])+1):
            ig = ig_['OBJ'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gSQR')
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y'+str(I)+','+str(int(v_['1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Y'+str(I)+','+str(int(v_['2']))])
            pbm.grelw = loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_['C1']
        posel = len(pbm.grelt[ig])
        pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['C1'])
        nlc = np.union1d(nlc,np.array([ig]))
        pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.007498464
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "SQR2-MN-4-1"
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
        pbm.efpar = arrset( pbm.efpar,0,np.sqrt(1.0e+0/6.2832e+0))
        return pbm

    @staticmethod
    def eY1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B = EV_[1]+EV_[2]*(1.0e+0-EV_[1])
        CI = pbm.elpar[iel_][0]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[0])
        P0V1 = -1.0e+0/(1.2e+1*EV_[0]**2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[0]**3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0**2)
        P1V1V1 = (2.0e+0*P0V1**2/P0-P0V1V1)/(P0**2)
        P2 = EV_[1]
        P3 = B**EV_[0]
        LOGB = np.log(B)
        P3V1 = P3*LOGB
        P3V2 = EV_[0]*(1.0e+0-EV_[2])*B**(EV_[0]-1.0e+0)
        P3V3 = EV_[0]*(1.0e+0-EV_[1])*B**(EV_[0]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*(1.0e+0-EV_[2])/B
        P3V1V3 = P3V3*LOGB+P3*(1.0e+0-EV_[1])/B
        P3V2V2 = EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[2])**2*B**(EV_[0]-1.0e+0)
        P3V2V3  = (
              -EV_[0]*B**(EV_[0]-1.0e+0)+EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[1])*(1.0e+0-EV_[2])*B**(EV_[0]-2.0e+0))
        P3V3V3 = EV_[0]*(EV_[0]-1.0e+0)*(1.0e+0-EV_[1])**2*B**(EV_[0]-2.0e+0)
        P4 = pbm.efpar[0]*np.sqrt(EV_[0])
        P4V1 = 5.0e-1*pbm.efpar[0]*np.sqrt(1.0e+0/EV_[0])
        P4V1V1 = -2.5e-1*pbm.efpar[0]*np.sqrt(1.0e+0/EV_[0]**3)
        C5 = CI**(-1.0e+0)
        P5 = C5*CI**EV_[0]
        P5V1 = P5*np.log(CI)
        P5V1V1 = P5V1*np.log(CI)
        P6 = np.exp(EV_[0]*(1.0e+0-CI*B))
        P6V1 = P6*(1.0e+0-CI*B)
        P6V2 = -P6*EV_[0]*CI*(1.0e+0-EV_[2])
        P6V3 = -P6*EV_[0]*CI*(1.0e+0-EV_[1])
        P6V1V1 = P6*(1.0e+0-CI*B)**2
        P6V1V2 = P6V2*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[2])
        P6V1V3 = P6V3*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_[1])
        P6V2V2 = -P6V2*EV_[0]*CI*(1.0e+0-EV_[2])
        P6V2V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[2])+P6*EV_[0]*CI
        P6V3V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[1])
        f_   = P1*P2*P3*P4*P5*P6
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[1] = P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[2] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[0,1]  = (
                      P1V1*(P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(P3V1*P4*P5*P6+(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[2,0] = H_[0,2]
                H_[1,1]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2+P3*P6V2+P3V2*P6)))
                H_[1,2]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3+P3V3*P6+P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[2,1] = H_[1,2]
                H_[2,2] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eY2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        B = EV_[1]+EV_[2]*(1.0e+0-EV_[1])
        CI = pbm.elpar[iel_][0]/7.658e+0
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_[0])
        P0V1 = -1.0e+0/(1.2e+1*EV_[0]**2)
        P0V1V1 = 2.0e+0/(1.2e+1*EV_[0]**3)
        P1 = 1.0e+0/P0
        P1V1 = -P0V1/(P0**2)
        P1V1V1 = (2.0e+0*P0V1**2/P0-P0V1V1)/(P0**2)
        P2 = 1.0e+0-EV_[1]
        P3 = (B/EV_[2])**EV_[0]
        LOGB = np.log(B/EV_[2])
        P3V1 = P3*LOGB
        P3V2 = EV_[0]*(-1.0e+0+1.0e+0/EV_[2])*(B/EV_[2])**(EV_[0]-1.0e+0)
        P3V3 = -EV_[0]*(EV_[1]/EV_[2]**2)*(B/EV_[2])**(EV_[0]-1.0e+0)
        P3V1V1 = P3V1*LOGB
        P3V1V2 = P3V2*LOGB+P3*EV_[2]*(-1.0e+0+1.0e+0/EV_[2])/B
        P3V1V3 = P3V3*LOGB-P3*EV_[1]/(B*EV_[2])
        P3V2V2  = (
              EV_[0]*(EV_[0]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[2])**2*(B/EV_[2])**(EV_[0]-2.0e+0))
        P3V2V3  = (
              EV_[0]*(-1.0e+0/EV_[2]**2)*(B/EV_[2])**(EV_[0]-1.0e+0)+EV_[0]*(EV_[0]-1.0e+0)*(-1.0e+0+1.0e+0/EV_[2])*(-EV_[1]/EV_[2]**2)*(B/EV_[2])**(EV_[0]-2.0e+0))
        P3V3V3 = (2.0e+0*EV_[0]*(EV_[1]/EV_[2]**3)*(B/EV_[2])**(EV_[0]-1.0e+0)+
             EV_[0]*(EV_[0]-1.0e+0)*(EV_[1]/EV_[2]**2)**2*(B/EV_[2])**(EV_[0]-2.0e+0))
        P4 = pbm.efpar[0]*np.sqrt(EV_[0])
        P4V1 = 5.0e-1*pbm.efpar[0]*np.sqrt(1.0e+0/EV_[0])
        P4V1V1 = -2.5e-1*pbm.efpar[0]*np.sqrt(1.0e+0/EV_[0]**3)
        C5 = CI**(-1.0e+0)
        P5 = C5*CI**EV_[0]
        P5V1 = P5*np.log(CI)
        P5V1V1 = P5V1*np.log(CI)
        P6 = np.exp(EV_[0]*(1.0e+0-CI*B/EV_[2]))
        P6V1 = P6*(1.0e+0-CI*B/EV_[2])
        P6V2 = -P6*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]
        P6V3 = P6*EV_[0]*CI*EV_[1]/EV_[2]**2
        P6V1V1 = P6*(1.0e+0-CI*B/EV_[2])**2
        P6V1V2 = P6V2*(1.0e+0-CI*B/EV_[2])-P6*CI*(-1.0e+0+1.0e+0/EV_[2])
        P6V1V3 = P6V3*(1.0e+0-CI*B/EV_[2])+P6*CI*EV_[1]/EV_[2]**2
        P6V2V2 = -P6V2*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]
        P6V2V3 = -P6V3*EV_[0]*CI*(1.0e+0-EV_[2])/EV_[2]+P6*EV_[0]*CI/EV_[2]**2
        P6V3V3 = (P6V3*EV_[0]*CI*EV_[1]/EV_[2]**2-2.0e+0*P6*EV_[0]*CI*EV_[1]/
             EV_[2]**3)
        f_   = P1*P2*P3*P4*P5*P6
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = (P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*P6+
                 P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1)
            g_[1] = -P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2
            g_[2] = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = (P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*P5*P6+
                     P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1))
                H_[0,1]  = (
                      P1V1*(-P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*(-P3V1*P4*P5*P6-(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2))))
                H_[1,0] = H_[0,1]
                H_[0,2]  = (
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3))))
                H_[2,0] = H_[0,2]
                H_[1,1]  = (
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2-P3*P6V2-P3V2*P6)))
                H_[1,2]  = (
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3-P3V3*P6-P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2))
                H_[2,1] = H_[1,2]
                H_[2,2] = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePROD(pbm,nargout,*args):

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
                H_[0,1] = 1.0e+0
                H_[1,0] = H_[0,1]
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
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

