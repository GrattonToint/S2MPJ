from s2mpjlib import *
class  CmRELOAD(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CmRELOAD
#    *********
# 
#    Source: Nuclear Reactor Core Reload Pattern Optimization
#    A.J. Quist et.al., draft paper, September 1997.
#    (2nd data set implemented here)
#    SIF input: S. Leyffer, November 1997
# 
#    classification = "C-CLOR2-MN-342-284"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CmRELOAD'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        v_['N'] = 14
        v_['M'] = 3
        v_['L'] = 4
        v_['T'] = 6
        v_['D11'] = 1
        v_['D12'] = 7
        v_['D21'] = 11
        v_['D22'] = 14
        v_['KFRESH'] = 1.25
        v_['FLIM'] = 1.8
        v_['KEFFuINI'] = 0.956145
        v_['ALPHA'] = 6E-6
        v_['CONSPW'] = 364.0
        v_['CYTIME'] = 350.0
        v_['TT'] = float(v_['T'])
        v_['T-1'] = -1.0+v_['TT']
        v_['DELTAT'] = v_['CYTIME']/v_['T-1']
        v_['ACC'] = v_['ALPHA']*v_['CONSPW']
        v_['ACC'] = v_['ACC']*v_['DELTAT']
        v_['-ACC'] = -1.0*v_['ACC']
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['V'+str(I)] = 1.0
        v_['V'+str(int(v_['D11']))] = 0.5
        v_['V'+str(int(v_['D12']))] = 0.5
        v_['V'+str(int(v_['D21']))] = 0.5
        v_['V'+str(int(v_['D22']))] = 0.5
        v_['NROW1'] = 2
        v_['ROWu'+str(int(v_['1']))+','+str(int(v_['1']))] = 1
        v_['ROWu'+str(int(v_['1']))+','+str(int(v_['2']))] = 2
        v_['G'+str(int(v_['1']))+','+str(int(v_['1']))] = 0.705
        v_['G'+str(int(v_['1']))+','+str(int(v_['2']))] = 0.25
        v_['NROW2'] = 5
        v_['ROWu'+str(int(v_['2']))+','+str(int(v_['1']))] = 1
        v_['ROWu'+str(int(v_['2']))+','+str(int(v_['2']))] = 2
        v_['ROWu'+str(int(v_['2']))+','+str(int(v_['3']))] = 3
        v_['ROWu'+str(int(v_['2']))+','+str(int(v_['4']))] = 7
        v_['ROWu'+str(int(v_['2']))+','+str(int(v_['5']))] = 8
        v_['G'+str(int(v_['2']))+','+str(int(v_['1']))] = 0.125
        v_['G'+str(int(v_['2']))+','+str(int(v_['2']))] = 0.625
        v_['G'+str(int(v_['2']))+','+str(int(v_['3']))] = 0.125
        v_['G'+str(int(v_['2']))+','+str(int(v_['7']))] = 0.08
        v_['G'+str(int(v_['2']))+','+str(int(v_['8']))] = 0.045
        v_['NROW3'] = 6
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['1']))] = 2
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['2']))] = 3
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['3']))] = 4
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['4']))] = 7
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['5']))] = 8
        v_['ROWu'+str(int(v_['3']))+','+str(int(v_['6']))] = 9
        v_['G'+str(int(v_['3']))+','+str(int(v_['2']))] = 0.125
        v_['G'+str(int(v_['3']))+','+str(int(v_['3']))] = 0.58
        v_['G'+str(int(v_['3']))+','+str(int(v_['4']))] = 0.125
        v_['G'+str(int(v_['3']))+','+str(int(v_['7']))] = 0.045
        v_['G'+str(int(v_['3']))+','+str(int(v_['8']))] = 0.08
        v_['G'+str(int(v_['3']))+','+str(int(v_['9']))] = 0.045
        v_['NROW4'] = 6
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['1']))] = 3
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['2']))] = 4
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['3']))] = 5
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['4']))] = 8
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['5']))] = 9
        v_['ROWu'+str(int(v_['4']))+','+str(int(v_['6']))] = 10
        v_['G'+str(int(v_['4']))+','+str(int(v_['3']))] = 0.125
        v_['G'+str(int(v_['4']))+','+str(int(v_['4']))] = 0.58
        v_['G'+str(int(v_['4']))+','+str(int(v_['5']))] = 0.125
        v_['G'+str(int(v_['4']))+','+str(int(v_['8']))] = 0.045
        v_['G'+str(int(v_['4']))+','+str(int(v_['9']))] = 0.08
        v_['G'+str(int(v_['4']))+','+str(int(v_['10']))] = 0.045
        v_['NROW5'] = 5
        v_['ROWu'+str(int(v_['5']))+','+str(int(v_['1']))] = 4
        v_['ROWu'+str(int(v_['5']))+','+str(int(v_['2']))] = 5
        v_['ROWu'+str(int(v_['5']))+','+str(int(v_['3']))] = 6
        v_['ROWu'+str(int(v_['5']))+','+str(int(v_['4']))] = 9
        v_['ROWu'+str(int(v_['5']))+','+str(int(v_['5']))] = 10
        v_['G'+str(int(v_['5']))+','+str(int(v_['4']))] = 0.125
        v_['G'+str(int(v_['5']))+','+str(int(v_['5']))] = 0.58
        v_['G'+str(int(v_['5']))+','+str(int(v_['6']))] = 0.125
        v_['G'+str(int(v_['5']))+','+str(int(v_['9']))] = 0.045
        v_['G'+str(int(v_['5']))+','+str(int(v_['10']))] = 0.08
        v_['NROW6'] = 3
        v_['ROWu'+str(int(v_['6']))+','+str(int(v_['1']))] = 5
        v_['ROWu'+str(int(v_['6']))+','+str(int(v_['2']))] = 6
        v_['ROWu'+str(int(v_['6']))+','+str(int(v_['3']))] = 10
        v_['G'+str(int(v_['6']))+','+str(int(v_['5']))] = 0.125
        v_['G'+str(int(v_['6']))+','+str(int(v_['6']))] = 0.58
        v_['G'+str(int(v_['6']))+','+str(int(v_['10']))] = 0.045
        v_['NROW7'] = 5
        v_['ROWu'+str(int(v_['7']))+','+str(int(v_['1']))] = 1
        v_['ROWu'+str(int(v_['7']))+','+str(int(v_['2']))] = 2
        v_['ROWu'+str(int(v_['7']))+','+str(int(v_['3']))] = 7
        v_['ROWu'+str(int(v_['7']))+','+str(int(v_['4']))] = 8
        v_['ROWu'+str(int(v_['7']))+','+str(int(v_['5']))] = 11
        v_['G'+str(int(v_['7']))+','+str(int(v_['1']))] = 0.045
        v_['G'+str(int(v_['7']))+','+str(int(v_['2']))] = 0.16
        v_['G'+str(int(v_['7']))+','+str(int(v_['7']))] = 0.5
        v_['G'+str(int(v_['7']))+','+str(int(v_['8']))] = 0.16
        v_['G'+str(int(v_['7']))+','+str(int(v_['11']))] = 0.045
        v_['NROW8'] = 8
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['1']))] = 2
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['2']))] = 3
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['3']))] = 4
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['4']))] = 7
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['5']))] = 8
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['6']))] = 9
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['7']))] = 11
        v_['ROWu'+str(int(v_['8']))+','+str(int(v_['8']))] = 12
        v_['G'+str(int(v_['8']))+','+str(int(v_['2']))] = 0.045
        v_['G'+str(int(v_['8']))+','+str(int(v_['3']))] = 0.08
        v_['G'+str(int(v_['8']))+','+str(int(v_['4']))] = 0.045
        v_['G'+str(int(v_['8']))+','+str(int(v_['7']))] = 0.08
        v_['G'+str(int(v_['8']))+','+str(int(v_['8']))] = 0.545
        v_['G'+str(int(v_['8']))+','+str(int(v_['9']))] = 0.08
        v_['G'+str(int(v_['8']))+','+str(int(v_['11']))] = 0.08
        v_['G'+str(int(v_['8']))+','+str(int(v_['12']))] = 0.045
        v_['NROW9'] = 9
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['1']))] = 3
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['2']))] = 4
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['3']))] = 5
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['4']))] = 8
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['5']))] = 9
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['6']))] = 10
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['7']))] = 11
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['8']))] = 12
        v_['ROWu'+str(int(v_['9']))+','+str(int(v_['9']))] = 13
        v_['G'+str(int(v_['9']))+','+str(int(v_['3']))] = 0.045
        v_['G'+str(int(v_['9']))+','+str(int(v_['4']))] = 0.08
        v_['G'+str(int(v_['9']))+','+str(int(v_['5']))] = 0.045
        v_['G'+str(int(v_['9']))+','+str(int(v_['8']))] = 0.08
        v_['G'+str(int(v_['9']))+','+str(int(v_['9']))] = 0.5
        v_['G'+str(int(v_['9']))+','+str(int(v_['10']))] = 0.08
        v_['G'+str(int(v_['9']))+','+str(int(v_['11']))] = 0.045
        v_['G'+str(int(v_['9']))+','+str(int(v_['12']))] = 0.08
        v_['G'+str(int(v_['9']))+','+str(int(v_['13']))] = 0.045
        v_['NROW10'] = 7
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['1']))] = 4
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['2']))] = 5
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['3']))] = 6
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['4']))] = 9
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['5']))] = 10
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['6']))] = 12
        v_['ROWu'+str(int(v_['10']))+','+str(int(v_['7']))] = 13
        v_['G'+str(int(v_['10']))+','+str(int(v_['4']))] = 0.045
        v_['G'+str(int(v_['10']))+','+str(int(v_['5']))] = 0.08
        v_['G'+str(int(v_['10']))+','+str(int(v_['6']))] = 0.045
        v_['G'+str(int(v_['10']))+','+str(int(v_['9']))] = 0.08
        v_['G'+str(int(v_['10']))+','+str(int(v_['10']))] = 0.5
        v_['G'+str(int(v_['10']))+','+str(int(v_['12']))] = 0.045
        v_['G'+str(int(v_['10']))+','+str(int(v_['13']))] = 0.08
        v_['NROW11'] = 6
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['1']))] = 7
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['2']))] = 8
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['3']))] = 9
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['4']))] = 11
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['5']))] = 12
        v_['ROWu'+str(int(v_['11']))+','+str(int(v_['6']))] = 14
        v_['G'+str(int(v_['11']))+','+str(int(v_['7']))] = 0.045
        v_['G'+str(int(v_['11']))+','+str(int(v_['8']))] = 0.125
        v_['G'+str(int(v_['11']))+','+str(int(v_['9']))] = 0.045
        v_['G'+str(int(v_['11']))+','+str(int(v_['11']))] = 0.5
        v_['G'+str(int(v_['11']))+','+str(int(v_['12']))] = 0.16
        v_['G'+str(int(v_['11']))+','+str(int(v_['14']))] = 0.045
        v_['NROW12'] = 7
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['1']))] = 8
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['2']))] = 9
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['3']))] = 10
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['4']))] = 11
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['5']))] = 12
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['6']))] = 13
        v_['ROWu'+str(int(v_['12']))+','+str(int(v_['7']))] = 14
        v_['G'+str(int(v_['12']))+','+str(int(v_['8']))] = 0.045
        v_['G'+str(int(v_['12']))+','+str(int(v_['9']))] = 0.08
        v_['G'+str(int(v_['12']))+','+str(int(v_['10']))] = 0.045
        v_['G'+str(int(v_['12']))+','+str(int(v_['11']))] = 0.08
        v_['G'+str(int(v_['12']))+','+str(int(v_['12']))] = 0.545
        v_['G'+str(int(v_['12']))+','+str(int(v_['13']))] = 0.08
        v_['G'+str(int(v_['12']))+','+str(int(v_['14']))] = 0.08
        v_['NROW13'] = 5
        v_['ROWu'+str(int(v_['13']))+','+str(int(v_['1']))] = 9
        v_['ROWu'+str(int(v_['13']))+','+str(int(v_['2']))] = 10
        v_['ROWu'+str(int(v_['13']))+','+str(int(v_['3']))] = 12
        v_['ROWu'+str(int(v_['13']))+','+str(int(v_['4']))] = 13
        v_['ROWu'+str(int(v_['13']))+','+str(int(v_['5']))] = 14
        v_['G'+str(int(v_['13']))+','+str(int(v_['9']))] = 0.045
        v_['G'+str(int(v_['13']))+','+str(int(v_['10']))] = 0.08
        v_['G'+str(int(v_['13']))+','+str(int(v_['12']))] = 0.08
        v_['G'+str(int(v_['13']))+','+str(int(v_['13']))] = 0.5
        v_['G'+str(int(v_['13']))+','+str(int(v_['14']))] = 0.045
        v_['NROW14'] = 3
        v_['ROWu'+str(int(v_['14']))+','+str(int(v_['1']))] = 11
        v_['ROWu'+str(int(v_['14']))+','+str(int(v_['2']))] = 12
        v_['ROWu'+str(int(v_['14']))+','+str(int(v_['3']))] = 14
        v_['G'+str(int(v_['14']))+','+str(int(v_['11']))] = 0.045
        v_['G'+str(int(v_['14']))+','+str(int(v_['12']))] = 0.125
        v_['G'+str(int(v_['14']))+','+str(int(v_['14']))] = 0.5
        v_['T-1'] = -1+v_['T']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(K)+','+str(J),ix_)
                    self.xnames=arrset(self.xnames,iv,'X'+str(I)+','+str(K)+','+str(J))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('KINF'+str(I)+','+str(S),ix_)
                self.xnames=arrset(self.xnames,iv,'KINF'+str(I)+','+str(S))
                [iv,ix_,_] = s2mpj_ii('PHI'+str(I)+','+str(S),ix_)
                self.xnames=arrset(self.xnames,iv,'PHI'+str(I)+','+str(S))
            [iv,ix_,_] = s2mpj_ii('KEFF'+str(S),ix_)
            self.xnames=arrset(self.xnames,iv,'KEFF'+str(S))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['KEFF'+str(int(v_['T']))]])
        valA = np.append(valA,float(-1.0))
        for K in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for I in range(int(v_['1']),int(v_['N'])+1):
                    [ig,ig_,_] = s2mpj_ii('SUMI'+str(K)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'SUMI'+str(K)+','+str(J))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(I)+','+str(K)+','+str(J)]])
                    valA = np.append(valA,float(v_['V'+str(I)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    [ig,ig_,_] = s2mpj_ii('SUMLM'+str(I),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'SUMLM'+str(I))
                    irA  = np.append(irA,[ig])
                    icA  = np.append(icA,[ix_['X'+str(I)+','+str(K)+','+str(J)]])
                    valA = np.append(valA,float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('PLAC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PLAC'+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['KINF'+str(I)+','+str(int(v_['1']))]])
            valA = np.append(valA,float(-1.0))
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('PLAC'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'PLAC'+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(I)+','+str(int(v_['1']))+','+str(J)]])
                valA = np.append(valA,float(v_['KFRESH']))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                [ig,ig_,_] = s2mpj_ii('KERN'+str(I)+','+str(S),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'KERN'+str(I)+','+str(S))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T-1'])+1):
                v_['R'] = 1+S
                [ig,ig_,_] = s2mpj_ii('KINFF'+str(I)+','+str(S),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'KINFF'+str(I)+','+str(S))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['KINF'+str(I)+','+str(int(v_['R']))]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['KINF'+str(I)+','+str(S)]])
                valA = np.append(valA,float(1.0))
        for S in range(int(v_['1']),int(v_['T'])+1):
            [ig,ig_,_] = s2mpj_ii('CPOW'+str(S),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'CPOW'+str(S))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                [ig,ig_,_] = s2mpj_ii('PEAK'+str(I)+','+str(S),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'PEAK'+str(I)+','+str(S))
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
        for K in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                self.gconst = arrset(self.gconst,ig_['SUMI'+str(K)+','+str(J)],float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            self.gconst = arrset(self.gconst,ig_['SUMLM'+str(I)],float(1.0))
        for S in range(int(v_['1']),int(v_['T'])+1):
            self.gconst = arrset(self.gconst,ig_['CPOW'+str(S)],float(1.0))
        v_['TEMP'] = 0.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['TEMP'] = v_['TEMP']+v_['V'+str(I)]
        v_['TEMP'] = v_['FLIM']/v_['TEMP']
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['PEAK'+str(I)+','+str(S)],float(v_['TEMP'])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    self.xupper[ix_['X'+str(I)+','+str(K)+','+str(J)]] = 1.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                self.xupper[ix_['KINF'+str(I)+','+str(S)]] = v_['KFRESH']
        v_['LOuKEFF'] = 0.9
        v_['UPuKEFF'] = 1.5
        v_['TEMP'] = float(v_['T'])
        v_['TEMP'] = -0.015*v_['TEMP']
        v_['LOuKEFF'] = v_['LOuKEFF']+v_['TEMP']
        v_['UPuKEFF'] = v_['UPuKEFF']+v_['TEMP']
        for S in range(int(v_['1']),int(v_['T'])+1):
            self.xlower[ix_['KEFF'+str(S)]] = v_['LOuKEFF']
            self.xupper[ix_['KEFF'+str(S)]] = v_['UPuKEFF']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        v_['R14'] = 14.0
        v_['TEMP'] = v_['KFRESH']/v_['R14']
        for S in range(int(v_['1']),int(v_['T'])+1):
            if('KEFF'+str(S) in ix_):
                self.x0[ix_['KEFF'+str(S)]] = float(v_['KEFFuINI'])
            else:
                self.y0  = (
                      arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['KEFF'+str(S)]),float(v_['KEFFuINI'])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                if('KINF'+str(I)+','+str(S) in ix_):
                    self.x0[ix_['KINF'+str(I)+','+str(S)]] = float(v_['KFRESH'])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['KINF'+str(I)+','+str(S)]),float(v_['KFRESH'])))
                if('PHI'+str(I)+','+str(S) in ix_):
                    self.x0[ix_['PHI'+str(I)+','+str(S)]] = float(v_['TEMP'])
                else:
                    self.y0  = (
                          arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['PHI'+str(I)+','+str(S)]),float(v_['TEMP'])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    if('X'+str(I)+','+str(K)+','+str(J) in ix_):
                        self.x0[ix_['X'+str(I)+','+str(K)+','+str(J)]] = float(0.5)
                    else:
                        self.y0  = (
                              arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['X'+str(I)+','+str(K)+','+str(J)]),float(0.5)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        [it,iet_,_] = s2mpj_ii( 'en3PROD', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        v_['K'] = 2
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Au'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype,ie,iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 3
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Bu'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype,ie,iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        v_['K'] = 4
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Cu'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    self.elftype = arrset(self.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype,ie,iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                    posev = np.where(elftv[ielftype[ie]]=='V3')[0]
                    self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ename = 'KTP'+str(I)+','+str(S)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype,ie,iet_["en2PROD"])
                vname = 'KEFF'+str(S)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'PHI'+str(I)+','+str(S)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ename = 'P'+str(I)+','+str(S)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype,ie,iet_["en2PROD"])
                vname = 'KINF'+str(I)+','+str(S)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'PHI'+str(I)+','+str(S)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
                posev = np.where(elftv[ielftype[ie]]=='V2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['PLAC'+str(I)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['Au'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,float(v_['V'+str(II)]))
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['Bu'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,float(v_['V'+str(II)]))
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['Cu'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw = loaset(self.grelw,ig,posel,float(v_['V'+str(II)]))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['KERN'+str(I)+','+str(S)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['KTP'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
                v_['TEMP'] = v_['NROW'+str(I)]
                v_['NuROW'] = int(np.fix(v_['TEMP']))
                for II in range(int(v_['1']),int(v_['NuROW'])+1):
                    v_['TEMP'] = v_['ROWu'+str(I)+','+str(II)]
                    v_['III'] = int(np.fix(v_['TEMP']))
                    ig = ig_['KERN'+str(I)+','+str(S)]
                    posel = len(self.grelt[ig])
                    self.grelt  = (
                          loaset(self.grelt,ig,posel,ie_['P'+str(int(v_['III']))+','+str(S)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    self.grelw  = (
                          loaset(self.grelw,ig,posel,float(v_['G'+str(I)+','+str(int(v_['III']))])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T-1'])+1):
                ig = ig_['KINFF'+str(I)+','+str(S)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['-ACC']))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['CPOW'+str(S)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['V'+str(I)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ig = ig_['PEAK'+str(I)+','+str(S)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.cupper[np.arange(self.nle)] = np.zeros((self.nle,1))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CLOR2-MN-342-284"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(self, nargout,*args):

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
    def en3PROD(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

