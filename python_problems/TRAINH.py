from s2mpjlib import *
class  TRAINH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TRAINH
#    *********
# 
#    The problem is to minimize the energy spent to move a train 
#    from the beginning of a track to its end in a given time.  The train
#    is slowed down by some drag (assumed to be quadratic in the the velocity).
#    The track follows the slope of a hill.  The track geometry is given
#    by the equation
# 
#           1                    1  ns-1                          x - z_i
#    g(x) = - ( a_1 + a_{ns} ) + -- SUM ( s_{i+1} - s_i ) arctan( ------- )
#           2                    pi  1                              eps
# 
#    where the z_i are the breakpoints between sections of the track, and where 
#    the s_i are the "slopes" on these sections (eps is a regularization
#    parameter). Here we have a track of the overall shape
# 
#                      ______
#                     /      \      z0 = 0, z1 = 2, z2 = 4, z3 = 6
#                    /        \     s1 = 2, s2 = 0, s3 = -2
#                   /          \    eps = 0.05
# 
#    The control variables are the acceleration force (UA) and the braking
#    force (UB) applied on the train.
# 
#    Source: adapted from
#    J. Kautsky and N. K. Nichols,
#    "OTEP-2: Optimal Train Energy Programme, mark 2",
#    Numerical Analysis Report NA/4/83,
#    Department of Mathematics, University of Reading, 1983.
# 
#    SIF input: N. Nichols and Ph. Toint, April 1993
# 
#    classification = "C-CQOR2-MN-V-V"
# 
#    Number of discretized points in the interval
# 
#           Alternative values for the SIF file parameters:
# IE N                   11             $-PARAMETER n=48, m=22
# IE N                   51             $-PARAMETER n=208, m=102
# IE N                   101            $-PARAMETER n=408, m=202  original value
# IE N                   201            $-PARAMETER n=808, m=402
# IE N                   501            $-PARAMETER n=2008, m=1002 
# IE N                   1001           $-PARAMETER n=4008, m=2002
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TRAINH'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(11);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   5001           $-PARAMETER n=20008, m=10002
        if nargin<2:
            v_['TIME'] = float(4.8);  #  SIF file default value
        else:
            v_['TIME'] = float(args[1])
        if nargin<3:
            v_['LENGTH'] = float(6.0);  #  SIF file default value
        else:
            v_['LENGTH'] = float(args[2])
        if nargin<4:
            v_['NS'] = int(3);  #  SIF file default value
        else:
            v_['NS'] = int(args[3])
        if nargin<5:
            v_['Z1'] = float(2.0);  #  SIF file default value
        else:
            v_['Z1'] = float(args[4])
        if nargin<6:
            v_['Z2'] = float(4.0);  #  SIF file default value
        else:
            v_['Z2'] = float(args[5])
        if nargin<7:
            v_['S1'] = float(2.0);  #  SIF file default value
        else:
            v_['S1'] = float(args[6])
        if nargin<8:
            v_['S2'] = float(0.0);  #  SIF file default value
        else:
            v_['S2'] = float(args[7])
        if nargin<9:
            v_['S3'] = float(-2.0);  #  SIF file default value
        else:
            v_['S3'] = float(args[8])
        v_['N-1'] = -1+v_['N']
        v_['RN'] = float(v_['N'])
        v_['H'] = v_['TIME']/v_['RN']
        v_['H/2'] = 0.5*v_['H']
        v_['-H'] = -1.0*v_['H']
        v_['-H/2'] = -1.0*v_['H/2']
        v_['UAMAX'] = 10.0
        v_['UBMIN'] = -2.0
        v_['VMAX'] = 10.0
        v_['A'] = 0.3
        v_['B'] = 0.14
        v_['C'] = 0.16
        v_['EPS'] = 0.05
        v_['0'] = 0
        v_['1'] = 1
        v_['PI'] = 3.1415926535
        v_['NS-1'] = -1+v_['NS']
        v_['BH/2'] = v_['B']*v_['H/2']
        v_['1+BH/2'] = 1.0+v_['BH/2']
        v_['BH/2-1'] = -1.0+v_['BH/2']
        v_['-AH'] = v_['A']*v_['-H']
        v_['LENGTH/N'] = v_['LENGTH']/v_['RN']
        v_['CH/2'] = v_['C']*v_['H/2']
        v_['SUMS'] = v_['S'+str(int(v_['1']))]+v_['S'+str(int(v_['NS']))]
        v_['-AVS'] = -0.5*v_['SUMS']
        v_['-AVSH'] = v_['-AVS']*v_['H']
        v_['CNST'] = v_['-AH']+v_['-AVSH']
        v_['H/2PI'] = v_['H/2']/v_['PI']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('V'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('UA'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'UA'+str(I))
        for I in range(int(v_['0']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('UB'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'UB'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('ENERGY',ig_)
        gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('XEQ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'XEQ'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['-H/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(I)]])
            valA = np.append(valA,float(v_['-H/2']))
            [ig,ig_,_] = s2mpj_ii('VEQ'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'VEQ'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['1+BH/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(I)]])
            valA = np.append(valA,float(v_['BH/2-1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UA'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['-H/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UA'+str(I)]])
            valA = np.append(valA,float(v_['-H/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UB'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['-H/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['UB'+str(I)]])
            valA = np.append(valA,float(v_['-H/2']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        legrps = np.where(gtype=='<=')[0]
        eqgrps = np.where(gtype=='==')[0]
        gegrps = np.where(gtype=='>=')[0]
        self.nle = len(legrps)
        self.neq = len(eqgrps)
        self.nge = len(gegrps)
        self.m   = self.nle+self.neq+self.nge
        self.congrps = np.concatenate((legrps,eqgrps,gegrps))
        self.cnames = cnames[self.congrps]
        self.nob = ngrp-self.m
        self.objgrps = np.where(gtype=='<>')[0]
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            self.gconst = arrset(self.gconst,ig_['VEQ'+str(I)],float(v_['CNST']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['X'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['V'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['V'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['UA'+str(int(v_['0']))]] = v_['UAMAX']
        self.xupper[ix_['UA'+str(int(v_['0']))]] = v_['UAMAX']
        self.xlower[ix_['UB'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['UB'+str(int(v_['0']))]] = 0.0
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            self.xlower[ix_['X'+str(I)]] = -float('Inf')
            self.xupper[ix_['X'+str(I)]] = +float('Inf')
            self.xlower[ix_['V'+str(I)]] = -float('Inf')
            self.xupper[ix_['V'+str(I)]] = +float('Inf')
            self.xlower[ix_['UA'+str(I)]] = 0.0
            self.xupper[ix_['UA'+str(I)]] = v_['UAMAX']
            self.xlower[ix_['UB'+str(I)]] = v_['UBMIN']
            self.xupper[ix_['UB'+str(I)]] = 0.0
        self.xlower[ix_['X'+str(int(v_['N']))]] = v_['LENGTH']
        self.xupper[ix_['X'+str(int(v_['N']))]] = v_['LENGTH']
        self.xlower[ix_['V'+str(int(v_['N']))]] = 0.0
        self.xupper[ix_['V'+str(int(v_['N']))]] = 0.0
        self.xlower[ix_['UA'+str(int(v_['N']))]] = 0.0
        self.xupper[ix_['UA'+str(int(v_['N']))]] = 0.0
        self.xlower[ix_['UB'+str(int(v_['N']))]] = v_['UBMIN']
        self.xupper[ix_['UB'+str(int(v_['N']))]] = v_['UBMIN']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['X'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['V'+str(int(v_['0']))]] = float(0.0)
        self.x0[ix_['UA'+str(int(v_['0']))]] = float(v_['UAMAX'])
        self.x0[ix_['UB'+str(int(v_['0']))]] = float(0.0)
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            v_['RI'] = float(I)
            v_['PI'] = v_['LENGTH/N']*v_['RI']
            self.x0[ix_['X'+str(I)]] = float(v_['PI'])
            self.x0[ix_['V'+str(I)]] = float(v_['LENGTH/N'])
            self.x0[ix_['UA'+str(I)]] = float(0.0)
            self.x0[ix_['UB'+str(I)]] = float(0.0)
        self.x0[ix_['X'+str(int(v_['N']))]] = float(v_['LENGTH'])
        self.x0[ix_['V'+str(int(v_['N']))]] = float(0.0)
        self.x0[ix_['UA'+str(int(v_['N']))]] = float(0.0)
        self.x0[ix_['UB'+str(int(v_['N']))]] = float(v_['UBMIN'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'UU')
        elftv = loaset(elftv,it,1,'VV')
        [it,iet_,_] = s2mpj_ii( 'eSQ', iet_)
        elftv = loaset(elftv,it,0,'VVV')
        [it,iet_,_] = s2mpj_ii( 'eATAN', iet_)
        elftv = loaset(elftv,it,0,'XX')
        elftp = []
        elftp = loaset(elftp,it,0,'ZZ')
        elftp = loaset(elftp,it,1,'E')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['0']),int(v_['N'])+1):
            v_['I+1'] = 1+I
            ename = 'VISQ'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSQ')
            ielftype = arrset(ielftype,ie,iet_["eSQ"])
            vname = 'V'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='VVV')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['1']),int(v_['NS-1'])+1):
                ename = 'A'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eATAN')
                ielftype = arrset(ielftype,ie,iet_["eATAN"])
                vname = 'X'+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='XX')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='ZZ')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['Z'+str(J)]))
                posep = np.where(elftp[ielftype[ie]]=='E')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['EPS']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ename = 'UV'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePROD')
            ielftype = arrset(ielftype,ie,iet_["ePROD"])
            vname = 'UA'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='UU')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'V'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='VV')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['N-1'])+1):
            v_['I+1'] = 1+I
            ig = ig_['VEQ'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['VISQ'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['CH/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['VISQ'+str(int(v_['I+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['CH/2']))
            for J in range(int(v_['1']),int(v_['NS-1'])+1):
                v_['J+1'] = 1+J
                v_['DS'] = v_['S'+str(int(v_['J+1']))]-v_['S'+str(J)]
                v_['WJ'] = v_['DS']*v_['H/2PI']
                ig = ig_['VEQ'+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['A'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WJ']))
                posel = len(self.grelt[ig])
                self.grelt  = (
                      loaset(self.grelt,ig,posel,ie_['A'+str(int(v_['I+1']))+','+str(J)]))
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['WJ']))
        for I in range(int(v_['1']),int(v_['N-1'])+1):
            ig = ig_['ENERGY']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['UV'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['H']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 0.0
# LO SOLUTION(11)        12.423025536
# LO SOLUTION(51)        12.295777964
# LO SOLUTION(101)       12.306399739
# LO SOLUTION(201)       12.309848614
# LO SOLUTION(1001)      12.307801327
# LO SOLUTION(5001)      12.221148056
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQOR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSQ(self, nargout,*args):

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
    def ePROD(self, nargout,*args):

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
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eATAN(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DX = EV_[0]-self.elpar[iel_][0]
        E2 = self.elpar[iel_][1]*self.elpar[iel_][1]
        f_   = np.arctan(DX/self.elpar[iel_][1])
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][1]/(E2+DX*DX)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -2.0*DX*self.elpar[iel_][1]/(E2+DX*DX)**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

