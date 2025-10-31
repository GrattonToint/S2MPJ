from s2mpjlib import *
class  CHEMRCTA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHEMRCTA
#    *********
# 
#    The tubular chemical reactor model problem by Poore, using a
#    finite difference approximation to the steady state solutions.
# 
#    Source: Problem 8, eqs (8.6)--(8.9) in
#    J.J. More',
#    "A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               minor correction by Ph. Shott, Jan 1995 and F Ruediger, Mar 1997.
# 
#    classification = "C-CNOR2-MN-V-V"
# 
#    The axial coordinate interval is [0,1]
# 
#    Number of discretized point for the interval [0,1].
#    The number of variables is 2N.
# 
#           Alternative values for the SIF file parameters:
# IE N                   5              $-PARAMETER n = 10
# IE N                   25             $-PARAMETER n = 50
# IE N                   50             $-PARAMETER n = 100
# IE N                   250            $-PARAMETER n = 500    original value
# IE N                   500            $-PARAMETER n = 1000
# IE N                   2500           $-PARAMETER n = 5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CHEMRCTA'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(5);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        if nargin<2:
            v_['PEM'] = float(1.0);  #  SIF file default value
        else:
            v_['PEM'] = float(args[1])
        if nargin<3:
            v_['PEH'] = float(5.0);  #  SIF file default value
        else:
            v_['PEH'] = float(args[2])
        if nargin<4:
            v_['D'] = float(0.135);  #  SIF file default value
        else:
            v_['D'] = float(args[3])
        if nargin<5:
            v_['B'] = float(0.5);  #  SIF file default value
        else:
            v_['B'] = float(args[4])
        if nargin<6:
            v_['BETA'] = float(2.0);  #  SIF file default value
        else:
            v_['BETA'] = float(args[5])
        if nargin<7:
            v_['GAMMA'] = float(25.0);  #  SIF file default value
        else:
            v_['GAMMA'] = float(args[6])
        v_['1'] = 1
        v_['2'] = 2
        v_['1.0'] = 1.0
        v_['N-1'] = -1+v_['N']
        v_['1/H'] = float(v_['N-1'])
        v_['-1/H'] = -1.0*v_['1/H']
        v_['H'] = v_['1.0']/v_['1/H']
        v_['1/H2'] = v_['1/H']*v_['1/H']
        v_['-D'] = -1.0*v_['D']
        v_['1/PEM'] = v_['1.0']/v_['PEM']
        v_['1/H2PEM'] = v_['1/PEM']*v_['1/H2']
        v_['-1/H2PM'] = -1.0*v_['1/H2PEM']
        v_['HPEM'] = v_['PEM']*v_['H']
        v_['-HPEM'] = -1.0*v_['HPEM']
        v_['-2/H2PM'] = v_['-1/H2PM']+v_['-1/H2PM']
        v_['CU1'] = 1.0*v_['-HPEM']
        v_['CUI-1'] = v_['1/H2PEM']+v_['1/H']
        v_['CUI'] = v_['-2/H2PM']+v_['-1/H']
        v_['BD'] = v_['B']*v_['D']
        v_['-BETA'] = -1.0*v_['BETA']
        v_['1/PEH'] = v_['1.0']/v_['PEH']
        v_['1/H2PEH'] = v_['1/PEH']*v_['1/H2']
        v_['-1/H2PH'] = -1.0*v_['1/H2PEH']
        v_['HPEH'] = v_['PEH']*v_['H']
        v_['-HPEH'] = -1.0*v_['HPEH']
        v_['-2/H2PH'] = v_['-1/H2PH']+v_['-1/H2PH']
        v_['CT1'] = 1.0*v_['-HPEH']
        v_['CTI-1'] = v_['1/H2PEH']+v_['1/H']
        v_['CTI'] = v_['-2/H2PH']+v_['-1/H']
        v_['CTI'] = v_['CTI']+v_['-BETA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('T'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'T'+str(I))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('U'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'U'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['2']))]])
        valA = np.append(valA,float(v_['CU1']))
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T'+str(int(v_['1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['1'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['1'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T'+str(int(v_['2']))]])
        valA = np.append(valA,float(v_['CT1']))
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            v_['I-1'] = -1+I
            v_['I+1'] = 1+I
            [ig,ig_,_] = s2mpj_ii('GU'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'GU'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['CUI-1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(I)]])
            valA = np.append(valA,float(v_['CUI']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['U'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['1/H2PEM']))
            [ig,ig_,_] = s2mpj_ii('GT'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'GT'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(v_['BETA']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(int(v_['I-1']))]])
            valA = np.append(valA,float(v_['CTI-1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)]])
            valA = np.append(valA,float(v_['CTI']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(int(v_['I+1']))]])
            valA = np.append(valA,float(v_['1/H2PEH']))
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N-1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('GU'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GU'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['U'+str(int(v_['N']))]])
        valA = np.append(valA,float(1.0))
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T'+str(int(v_['N-1']))]])
        valA = np.append(valA,float(-1.0))
        [ig,ig_,_] = s2mpj_ii('GT'+str(int(v_['N'])),ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'GT'+str(int(v_['N'])))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['T'+str(int(v_['N']))]])
        valA = np.append(valA,float(1.0))
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
        self.gconst  = (
              arrset(self.gconst,ig_['GU'+str(int(v_['1']))],float(v_['-HPEM'])))
        self.gconst  = (
              arrset(self.gconst,ig_['GT'+str(int(v_['1']))],float(v_['-HPEH'])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xlower[ix_['T'+str(I)]] = 0.0000001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eREAC', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'T')
        elftp = []
        elftp = loaset(elftp,it,0,'G')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            ename = 'EU'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eREAC')
                ielftype = arrset(ielftype,ie,iet_['eREAC'])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='G')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['GAMMA']))
            ename = 'ET'+str(I)
            [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
            if newelt:
                self.elftype = arrset(self.elftype,ie,'eREAC')
                ielftype = arrset(ielftype,ie,iet_['eREAC'])
            vname = 'U'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='U')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'T'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='T')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='G')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['GAMMA']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['2']),int(v_['N-1'])+1):
            ig = ig_['GU'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EU'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-D']))
            ig = ig_['GT'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['ET'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['BD']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
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
        self.pbclass   = "C-CNOR2-MN-V-V"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eREAC(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        DADT = self.elpar[iel_][0]/(EV_[1]*EV_[1])
        D2ADT2 = -2.0*DADT/EV_[1]
        EX = np.exp(self.elpar[iel_][0]-self.elpar[iel_][0]/EV_[1])
        UEX = EX*EV_[0]
        f_   = UEX
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EX
            g_[1] = UEX*DADT
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = EX*DADT
                H_[1,0] = H_[0,1]
                H_[1,1] = UEX*(DADT*DADT+D2ADT2)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

