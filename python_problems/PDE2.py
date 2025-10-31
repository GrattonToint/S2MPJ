from s2mpjlib import *
class  PDE2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PDE2
#    *********
# 
#    The pde_2, _20 & _200.mod AMPL models from Hans Mittelmann 
#    (mittelmann@asu.edu)
#    See: http://plato.asu.edu/ftp/barrier/
# 
#    SIF input: Nick Gould, April 25th 2012
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CLLR2-AN-V-V"
# 
#    the x-y discretization 
# 
#           Alternative values for the SIF file parameters:
# IE N                   3              $-PARAMETER
# IE N                   299            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PDE2'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(6);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
# IE N                   2999           $-PARAMETER     pde_2.mod value
# IE N                   2099           $-PARAMETER     pde_20.mod value
# IE N                   1299           $-PARAMETER     pde_200.mod value
        v_['0'] = 0
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['N1'] = 1+v_['N']
        v_['RN1'] = float(v_['N1'])
        v_['A'] = 0.01
        v_['G'] = 20.0
        v_['H'] = v_['ONE']/v_['RN1']
        v_['-H'] = -1.0*v_['H']
        v_['H2'] = v_['H']*v_['H']
        v_['GH2'] = v_['G']*v_['H2']
        v_['AH'] = v_['A']*v_['H']
        v_['SQRTAH'] = np.sqrt(v_['AH'])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['0']),int(v_['N1'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2mpj_ii('T'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'T'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+'] = 1+I
            v_['I-'] = -1+I
            for J in range(int(v_['1']),int(v_['N'])+1):
                v_['J+'] = 1+J
                v_['J-'] = -1+J
                [ig,ig_,_] = s2mpj_ii('P'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'P'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(4.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['J+']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['J-']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I+']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I-']))+','+str(J)]])
                valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(I)+','+str(J)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['H']))
                [ig,ig_,_] = s2mpj_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'B'+str(I)+','+str(J))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(I)+','+str(J)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(J)]])
                valA = np.append(valA,float(v_['H']))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['0'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['0'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['0'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['0'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['N1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['N1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['N1']))]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['N1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['N1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['N1']))]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['0']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['0']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['0']))+','+str(I)]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['0']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['0']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['0']))+','+str(I)]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['N1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['N1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['N1']))+','+str(I)]])
            valA = np.append(valA,float(v_['SQRTAH']))
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['N1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(I)+','+str(int(v_['0']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['N1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['X'+str(int(v_['N1']))+','+str(I)]])
            valA = np.append(valA,float(v_['SQRTAH']))
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['P'+str(I)+','+str(J)],float(v_['GH2'])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['RI'] = float(I)
            v_['IH'] = v_['RI']*v_['H']
            v_['IH-1'] = -1.0+v_['IH']
            v_['P'] = v_['RI']*v_['IH-1']
            v_['P'] = 5.0*v_['P']
            for J in range(int(v_['1']),int(v_['N'])+1):
                v_['RJ'] = float(J)
                v_['JH'] = v_['RJ']*v_['H']
                v_['JH-1'] = -1.0+v_['JH']
                v_['YD'] = v_['RJ']*v_['JH-1']
                v_['YD'] = v_['YD']*v_['P']
                v_['YD'] = 3.0+v_['YD']
                v_['YD'] = v_['YD']*v_['H']
                v_['-YD'] = -1.0*v_['YD']
                self.gconst = arrset(self.gconst,ig_['A'+str(I)+','+str(J)],float(v_['YD']))
                self.gconst  = (
                      arrset(self.gconst,ig_['B'+str(I)+','+str(J)],float(v_['-YD'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),0.0)
        self.xupper = np.full((self.n,1),3.5)
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.xupper[ix_['X'+str(I)+','+str(int(v_['0']))]] = 10.0
            self.xupper[ix_['X'+str(I)+','+str(int(v_['N1']))]] = 10.0
            self.xupper[ix_['X'+str(int(v_['0']))+','+str(I)]] = 10.0
            self.xupper[ix_['X'+str(int(v_['N1']))+','+str(I)]] = 10.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons   = np.arange(len(self.congrps))
        self.pbclass   = "C-CLLR2-AN-V-V"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2
        self.conderlvl = [2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

