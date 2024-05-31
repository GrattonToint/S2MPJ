from s2xlib import *
class  MANCINO(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MANCINO
#    *********
# 
#    Mancino's function with variable dimension.
# 
#    Source:
#    E. Spedicato,
#    "Computational experience with quasi-Newton algorithms for
#    minimization problems of moderate size",
#    Report N-175, CISE, Milano, 1975.
# 
#    See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January, 1995.
# 
#    classification = "SUR2-AN-V-0"
# 
#    The definitions
#      s_{i,j} = \sin \log v_{i,j}   and s_{i,j} = \cos \log v_{i,j}
#    have been used.  It seems that the additional exponent ALPHA
#    in Buckley is a typo.
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MANCINO'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MANCINO'
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
# IE N                   20             $-PARAMETER
# IE N                   30             $-PARAMETER Schittkowski #391
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
        if nargin<2:
            v_['ALPHA'] = int(5);  #  SIF file default value
        else:
            v_['ALPHA'] = int(args[1])
        if nargin<3:
            v_['BETA'] = float(14.0);  #  SIF file default value
        else:
            v_['BETA'] = float(args[2])
        if nargin<4:
            v_['GAMMA'] = int(3);  #  SIF file default value
        else:
            v_['GAMMA'] = int(args[3])
        v_['RALPHA'] = float(v_['ALPHA'])
        v_['RN'] = float(v_['N'])
        v_['N-1'] = -1+v_['N']
        v_['RN-1'] = float(v_['N-1'])
        v_['N-1SQ'] = v_['RN-1']*v_['RN-1']
        v_['BETAN'] = v_['BETA']*v_['RN']
        v_['BETAN2'] = v_['BETAN']*v_['BETAN']
        v_['AL+1'] = 1.0+v_['RALPHA']
        v_['A1SQ'] = v_['AL+1']*v_['AL+1']
        v_['F0'] = v_['A1SQ']*v_['N-1SQ']
        v_['F1'] = -1.0*v_['F0']
        v_['F2'] = v_['BETAN2']+v_['F1']
        v_['F3'] = 1.0/v_['F2']
        v_['F4'] = v_['BETAN']*v_['F3']
        v_['A'] = -1.0*v_['F4']
        v_['-N/2'] = -0.5*v_['RN']
        v_['1'] = 1
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
            [ig,ig_,_] = s2x_ii('G'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
            iv = ix_['X'+str(I)]
            pbm.A[ig,iv] = float(v_['BETAN'])+pbm.A[ig,iv]
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['I-N/2'] = v_['RI']+v_['-N/2']
            v_['CI'] = 1.0
            for J in range(int(v_['1']),int(v_['GAMMA'])+1):
                v_['CI'] = v_['CI']*v_['I-N/2']
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)],float(v_['CI']))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            v_['RI'] = float(I)
            v_['H'] = 0.0
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RJ'] = float(J)
                v_['1/J'] = 1.0/v_['RJ']
                v_['I/J'] = v_['RI']*v_['1/J']
                v_['SQI/J'] = np.sqrt(v_['I/J'])
                v_['LIJ'] = np.log(v_['SQI/J'])
                v_['SIJ'] = np.sin(v_['LIJ'])
                v_['CIJ'] = np.cos(v_['LIJ'])
                v_['SA'] = 1.0
                v_['CA'] = 1.0
                for K in range(int(v_['1']),int(v_['ALPHA'])+1):
                    v_['SA'] = v_['SA']*v_['SIJ']
                    v_['CA'] = v_['CA']*v_['CIJ']
                v_['SCA'] = v_['SA']+v_['CA']
                v_['HIJ'] = v_['SQI/J']*v_['SCA']
                v_['H'] = v_['H']+v_['HIJ']
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['1/J'] = 1.0/v_['RJ']
                v_['I/J'] = v_['RI']*v_['1/J']
                v_['SQI/J'] = np.sqrt(v_['I/J'])
                v_['LIJ'] = np.log(v_['SQI/J'])
                v_['SIJ'] = np.sin(v_['LIJ'])
                v_['CIJ'] = np.cos(v_['LIJ'])
                v_['SA'] = 1.0
                v_['CA'] = 1.0
                for K in range(int(v_['1']),int(v_['ALPHA'])+1):
                    v_['SA'] = v_['SA']*v_['SIJ']
                    v_['CA'] = v_['CA']*v_['CIJ']
                v_['SCA'] = v_['SA']+v_['CA']
                v_['HIJ'] = v_['SQI/J']*v_['SCA']
                v_['H'] = v_['H']+v_['HIJ']
            v_['I-N/2'] = v_['RI']+v_['-N/2']
            v_['CI'] = 1.0
            for J in range(int(v_['1']),int(v_['GAMMA'])+1):
                v_['CI'] = v_['CI']*v_['I-N/2']
            v_['TMP'] = v_['H']+v_['CI']
            v_['XI0'] = v_['TMP']*v_['A']
            pb.x0[ix_['X'+str(I)]] = float(v_['XI0'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eMANC', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftp = []
        elftp = loaset(elftp,it,0,'II')
        elftp = loaset(elftp,it,1,'JJ')
        elftp = loaset(elftp,it,2,'AL')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                v_['RJ'] = float(J)
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eMANC')
                ielftype = arrset(ielftype, ie, iet_["eMANC"])
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='II')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RI']))
                posep = find(elftp[ielftype[ie]],lambda x:x=='JJ')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RJ']))
                posep = find(elftp[ielftype[ie]],lambda x:x=='AL')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RALPHA']))
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eMANC')
                ielftype = arrset(ielftype, ie, iet_["eMANC"])
                vname = 'X'+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                posep = find(elftp[ielftype[ie]],lambda x:x=='II')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RI']))
                posep = find(elftp[ielftype[ie]],lambda x:x=='JJ')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RJ']))
                posep = find(elftp[ielftype[ie]],lambda x:x=='AL')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['RALPHA']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I-1'] = -1+I
            ig = ig_['G'+str(I)]
            pbm.grftype = arrset(pbm.grftype,ig,'gL2')
            for J in range(int(v_['1']),int(v_['I-1'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            v_['I+1'] = 1+I
            for J in range(int(v_['I+1']),int(v_['N'])+1):
                ig = ig_['G'+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SUR2-AN-V-0"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eMANC(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        IAL = pbm.elpar[iel_][2]
        IA1 = IAL-1
        A2 = pbm.elpar[iel_][2]-2.0
        IA2 = IAL-2
        IA3 = IAL-3
        INVIJ = EV_[0]*EV_[0]+pbm.elpar[iel_][0]/pbm.elpar[iel_][1]
        VIJ = np.sqrt(INVIJ)
        V2 = VIJ*VIJ
        DVIJ = EV_[0]/VIJ
        LIJ = np.log(VIJ)
        SIJ = np.sin(LIJ)
        CIJ = np.cos(LIJ)
        DSDX = CIJ*DVIJ/VIJ
        DCDX = -SIJ*DVIJ/VIJ
        SUMAL = SIJ**IAL+CIJ**IAL
        DSUMAL = pbm.elpar[iel_][2]*(DSDX*SIJ**IA1+DCDX*CIJ**IA1)
        SCIJ = SIJ*CIJ
        DSCIJ = SIJ*DCDX+DSDX*CIJ
        SAL = SIJ**IA2-CIJ**IA2
        DSAL = A2*(DSDX*SIJ**IA3-DCDX*CIJ**IA3)
        B = SUMAL+pbm.elpar[iel_][2]*SCIJ*SAL
        DBDX = DSUMAL+pbm.elpar[iel_][2]*(DSCIJ*SAL+SCIJ*DSAL)
        f_   = VIJ*SUMAL
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]*B/VIJ
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = (B+EV_[0]*DBDX)/VIJ-B*EV_[0]*DVIJ/V2
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

