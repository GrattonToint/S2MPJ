from s2mpjlib import *
class  DRUGDIS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DRUGDIS
#    *********
# 
#    A control problem based on the kinetic model of Aarons and Rowland for
#    DRUG DISplacemnt, which simulates the interaction of the two drugs 
#    (warfarin and phenylnutazone) in a patient bloodstream.  
#    The state variable are the concentrations of unbound warfarin (w) and 
#    phenylbutazone (p).  The problem is to control the rate of injection (u) 
#    of the pain-killing phenylbutazone so that both drugs reach a specified 
#    steady-state in minimum time and the concentration of warfarin does not 
#    rise above a given toxicity level.
# 
#    The problem is discretized using the trapezoidal rule.  It is non-convex.
# 
#    The problem can be made harder by diminishing the value of the lower bound
#    on the final time TF (while maintaining it strictly positive).
# 
#    Source:
#    H. Maurer and M. Wiegand,
#    "Numerical solution of a drug displacement problem with bounded state
#    variables",
#    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
# 
#    SIF input: Ph. Toint, Nov 1993.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "LOR2-MN-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#           Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=  34, m= 20 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'DRUGDIS'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'DRUGDIS'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NI'] = int(10);  #  SIF file default value
        else:
            v_['NI'] = int(args[0])
# IE NI                  50             $-PARAMETER n= 154, m=100 
# IE NI                  100            $-PARAMETER n= 304, m=200  original value
# IE NI                  200            $-PARAMETER n= 604, m=400 
# IE NI                  500            $-PARAMETER n=1504, m=1000 
# IE NI                  1000           $-PARAMETER n=3004, m=2000 
# IE NI                  2000           $-PARAMETER n=6004, m=4000 
        if nargin<2:
            v_['TOXIC'] = float(0.026);  #  SIF file default value
        else:
            v_['TOXIC'] = float(args[1])
        if nargin<3:
            v_['WSS'] = float(0.02);  #  SIF file default value
        else:
            v_['WSS'] = float(args[2])
        if nargin<4:
            v_['UMAX'] = float(8.0);  #  SIF file default value
        else:
            v_['UMAX'] = float(args[3])
        if nargin<5:
            v_['PSTART'] = float(0.0);  #  SIF file default value
        else:
            v_['PSTART'] = float(args[4])
        if nargin<6:
            v_['PFINAL'] = float(2.0);  #  SIF file default value
        else:
            v_['PFINAL'] = float(args[5])
        v_['AVP'] = v_['PSTART']+v_['PFINAL']
        v_['AVP'] = 0.5*v_['AVP']
        v_['NI-1'] = -1+v_['NI']
        v_['RNI'] = float(v_['NI'])
        v_['-1/2NI'] = -0.5/v_['RNI']
        v_['0'] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('TF',ix_)
        pb.xnames=arrset(pb.xnames,iv,'TF')
        xscale = arrset(xscale,iv,200.0)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('W'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I))
            xscale = arrset(xscale,iv,0.02)
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('P'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'P'+str(I))
        for I in range(int(v_['0']),int(v_['NI'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('TFINAL',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['TF']
        pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        pbm.gscale = arrset(pbm.gscale,ig,float(100.0))
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('EW'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EW'+str(I))
            iv = ix_['W'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['W'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            pbm.gscale = arrset(pbm.gscale,ig,float(0.02))
            [ig,ig_,_] = s2mpj_ii('EP'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EP'+str(I))
            iv = ix_['P'+str(int(v_['I+1']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['P'+str(I)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['TF']] = 200.0
        for I in range(int(v_['0']),int(v_['NI'])+1):
            pb.xupper[ix_['W'+str(I)]] = v_['TOXIC']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            pb.xupper[ix_['U'+str(I)]] = v_['UMAX']
        pb.xlower[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        pb.xupper[ix_['W'+str(int(v_['0']))]] = v_['WSS']
        pb.xlower[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        pb.xupper[ix_['W'+str(int(v_['NI']))]] = v_['WSS']
        pb.xlower[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        pb.xupper[ix_['P'+str(int(v_['0']))]] = v_['PSTART']
        pb.xlower[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        pb.xupper[ix_['P'+str(int(v_['NI']))]] = v_['PFINAL']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        v_['DP'] = v_['PFINAL']-v_['PSTART']
        v_['DP/NI'] = v_['DP']/v_['RNI']
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            v_['RI'] = float(I)
            v_['IDP/NI'] = v_['RI']*v_['DP/NI']
            pb.x0[ix_['P'+str(I)]] = float(v_['IDP/NI'])
            pb.x0[ix_['W'+str(I)]] = float(v_['WSS'])
            pb.x0[ix_['U'+str(I)]] = float(v_['UMAX'])
        pb.x0[ix_['TF']] = float(240.0)
        pb.x0[ix_['W'+str(int(v_['NI']))]] = float(v_['WSS'])
        pb.x0[ix_['P'+str(int(v_['NI']))]] = float(v_['PFINAL'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEW', iet_)
        elftv = loaset(elftv,it,0,'T')
        elftv = loaset(elftv,it,1,'W')
        elftv = loaset(elftv,it,2,'P')
        elftv = loaset(elftv,it,3,'U')
        [it,iet_,_] = s2mpj_ii( 'eEP', iet_)
        elftv = loaset(elftv,it,0,'T')
        elftv = loaset(elftv,it,1,'W')
        elftv = loaset(elftv,it,2,'P')
        elftv = loaset(elftv,it,3,'U')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['NI'])+1):
            ename = 'WA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEW')
            ielftype = arrset(ielftype, ie, iet_["eEW"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='U')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'PA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eEP')
            ielftype = arrset(ielftype, ie, iet_["eEP"])
            vname = 'TF'
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='T')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'W'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'P'+str(I)
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='P')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            vname = 'U'+str(I)
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
        for I in range(int(v_['0']),int(v_['NI-1'])+1):
            ig = ig_['EW'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['WA'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/2NI']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['WA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/2NI']))
            ig = ig_['EP'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PA'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/2NI']))
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['PA'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-1/2NI']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 200.0
#    Solution
# LO SOLTN(10)           3.82432
# LO SOLTN(50)           4.19953
# LO SOLTN(100)          4.23934
# LO SOLTN(200)          4.25762
# LO SOLTN(500)
# LO SOLTN(1000)
# LO SOLTN(Maurer)       2.62637
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
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
        #%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        lxs = len(xscale);
        for j in np.arange(0,min(sA2,pb.n,len(xscale))):
            if not xscale[j] is None and xscale[j] != 0.0 and xscale[j] != 1.0:
                for i in find(pbm.A[:,j],lambda x:x!=0):
                      pbm.A[i,j] = pbm.A[i,j]/xscale[j]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-MN-V-V"
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
        pbm.efpar = arrset( pbm.efpar,0,46.4)
        pbm.efpar = arrset( pbm.efpar,1,0.02)
        pbm.efpar = arrset( pbm.efpar,2,0.2)
        pbm.efpar = arrset( pbm.efpar,3,232.0)
        pbm.efpar = arrset( pbm.efpar,4,pbm.efpar[0]*pbm.efpar[0])
        return pbm

    @staticmethod
    def eEW(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        D = 1.0+pbm.efpar[2]*(EV_[1]+EV_[2])
        DD = D*D
        DD1 = 2.0*pbm.efpar[2]*D
        DD2 = 2.0*pbm.efpar[2]*pbm.efpar[2]
        A = DD+pbm.efpar[3]+pbm.efpar[0]*EV_[1]
        AW = DD1+pbm.efpar[0]
        B = DD+pbm.efpar[3]+pbm.efpar[0]*EV_[2]
        BP = DD1+pbm.efpar[0]
        C = A*B-pbm.efpar[4]*EV_[1]*EV_[2]
        CW = AW*B+A*DD1-pbm.efpar[4]*EV_[2]
        CP = DD1*B+A*BP-pbm.efpar[4]*EV_[1]
        CWW = DD2*B+2.0*AW*DD1+A*DD2
        CWP = DD2*B+AW*BP+DD1*DD1+A*DD2-pbm.efpar[4]
        CPP = DD2*B+2.0*DD1*BP+A*DD2
        F = DD/C
        H = DD1-F*CW
        I = DD1-F*CP
        FW = H/C
        FP = I/C
        HW = DD2-CW*FW-F*CWW
        HP = DD2-CW*FW-F*CWP
        IP = DD2-CP*FP-F*CPP
        FWW = (HW-FW*CW)/C
        FWP = (HP-FW*CP)/C
        FPP = (IP-FP*CP)/C
        GU = pbm.efpar[0]*EV_[1]
        G = A*(pbm.efpar[1]-EV_[1])+GU*(EV_[3]-2.0*EV_[2])
        GW = AW*(pbm.efpar[1]-EV_[1])-A+pbm.efpar[0]*(EV_[3]-2.0*EV_[2])
        GP = DD1*(pbm.efpar[1]-EV_[1])-2.0*GU
        GPP = DD2*(pbm.efpar[1]-EV_[1])
        GWW = GPP-2.0*AW
        GWP = GPP-DD1-2.0*pbm.efpar[0]
        f_   = EV_[0]*F*G
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = F*G
            g_[1] = EV_[0]*(FW*G+F*GW)
            g_[2] = EV_[0]*(FP*G+F*GP)
            g_[3] = EV_[0]*F*GU
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = FW*G+F*GW
                H_[1,0] = H_[0,1]
                H_[0,2] = FP*G+F*GP
                H_[2,0] = H_[0,2]
                H_[0,3] = F*GU
                H_[3,0] = H_[0,3]
                H_[1,1] = EV_[0]*(FWW*G+2.0*FW*GW+F*GWW)
                H_[1,2] = EV_[0]*(FWP*G+FW*GP+FP*GW+F*GWP)
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*(FW*GU+F*pbm.efpar[0])
                H_[3,1] = H_[1,3]
                H_[2,2] = EV_[0]*(FPP*G+2.0*FP*GP+F*GPP)
                H_[2,3] = EV_[0]*FP*GU
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eEP(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        D = 1.0+pbm.efpar[2]*(EV_[1]+EV_[2])
        DD = D*D
        DD1 = 2.0*pbm.efpar[2]*D
        DD2 = 2.0*pbm.efpar[2]*pbm.efpar[2]
        A = DD+pbm.efpar[3]+pbm.efpar[0]*EV_[1]
        AW = DD1+pbm.efpar[0]
        B = DD+pbm.efpar[3]+pbm.efpar[0]*EV_[2]
        BP = DD1+pbm.efpar[0]
        C = A*B-pbm.efpar[4]*EV_[1]*EV_[2]
        CW = AW*B+A*DD1-pbm.efpar[4]*EV_[2]
        CP = DD1*B+A*BP-pbm.efpar[4]*EV_[1]
        CWW = DD2*B+2.0*AW*DD1+A*DD2
        CWP = DD2*B+AW*BP+DD1*DD1+A*DD2-pbm.efpar[4]
        CPP = DD2*B+2.0*DD1*BP+A*DD2
        F = DD/C
        H = DD1-F*CW
        I = DD1-F*CP
        FW = H/C
        FP = I/C
        HW = DD2-CW*FW-F*CWW
        HP = DD2-CW*FW-F*CWP
        IP = DD2-CP*FP-F*CPP
        FWW = (HW-FW*CW)/C
        FWP = (HP-FW*CP)/C
        FPP = (IP-FP*CP)/C
        G = B*(EV_[3]-2.0*EV_[2])+pbm.efpar[0]*EV_[2]*(pbm.efpar[1]-EV_[1])
        GW = DD1*(EV_[3]-2.0*EV_[2])-pbm.efpar[0]*EV_[2]
        GP = BP*(EV_[3]-2.0*EV_[2])-2.0*B+pbm.efpar[0]*(pbm.efpar[1]-EV_[1])
        GWW = DD2*(EV_[3]-2.0*EV_[2])
        GWP = GWW-2.0*DD1-pbm.efpar[0]
        GPP = GWW-4.0*BP
        f_   = EV_[0]*F*G
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = F*G
            g_[1] = EV_[0]*(FW*G+F*GW)
            g_[2] = EV_[0]*(FP*G+F*GP)
            g_[3] = EV_[0]*F*B
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,1] = FW*G+F*GW
                H_[1,0] = H_[0,1]
                H_[0,2] = FP*G+F*GP
                H_[2,0] = H_[0,2]
                H_[0,3] = F*B
                H_[3,0] = H_[0,3]
                H_[1,1] = EV_[0]*(FWW*G+2.0*FW*GW+F*GWW)
                H_[1,2] = EV_[0]*(FWP*G+FW*GP+FP*GW+F*GWP)
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[0]*(FW*B+F*DD1)
                H_[3,1] = H_[1,3]
                H_[2,2] = EV_[0]*(FPP*G+2.0*FP*GP+F*GPP)
                H_[2,3] = EV_[0]*(FP*B+F*BP)
                H_[3,2] = H_[2,3]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

