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
#    classification = "LOR2-MN-342-284"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'CmRELOAD'

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
        pb.xnames = np.array([])
        pb.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(K)+','+str(J),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(K)+','+str(J))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('KINF'+str(I)+','+str(S),ix_)
                pb.xnames=arrset(pb.xnames,iv,'KINF'+str(I)+','+str(S))
                [iv,ix_,_] = s2mpj_ii('PHI'+str(I)+','+str(S),ix_)
                pb.xnames=arrset(pb.xnames,iv,'PHI'+str(I)+','+str(S))
            [iv,ix_,_] = s2mpj_ii('KEFF'+str(S),ix_)
            pb.xnames=arrset(pb.xnames,iv,'KEFF'+str(S))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['KEFF'+str(int(v_['T']))]
        pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for K in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for I in range(int(v_['1']),int(v_['N'])+1):
                    [ig,ig_,_] = s2mpj_ii('SUMI'+str(K)+','+str(J),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'SUMI'+str(K)+','+str(J))
                    iv = ix_['X'+str(I)+','+str(K)+','+str(J)]
                    pbm.A[ig,iv] = float(v_['V'+str(I)])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    [ig,ig_,_] = s2mpj_ii('SUMLM'+str(I),ig_)
                    gtype = arrset(gtype,ig,'==')
                    cnames = arrset(cnames,ig,'SUMLM'+str(I))
                    iv = ix_['X'+str(I)+','+str(K)+','+str(J)]
                    pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('PLAC'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'PLAC'+str(I))
            iv = ix_['KINF'+str(I)+','+str(int(v_['1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for J in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('PLAC'+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'PLAC'+str(I))
                iv = ix_['X'+str(I)+','+str(int(v_['1']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['KFRESH'])+pbm.A[ig,iv]
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
                iv = ix_['KINF'+str(I)+','+str(int(v_['R']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['KINF'+str(I)+','+str(S)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        for K in range(int(v_['1']),int(v_['L'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['SUMI'+str(K)+','+str(J)],float(1.0))
        for I in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['SUMLM'+str(I)],float(1.0))
        for S in range(int(v_['1']),int(v_['T'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['CPOW'+str(S)],float(1.0))
        v_['TEMP'] = 0.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['TEMP'] = v_['TEMP']+v_['V'+str(I)]
        v_['TEMP'] = v_['FLIM']/v_['TEMP']
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['PEAK'+str(I)+','+str(S)],float(v_['TEMP'])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    pb.xupper[ix_['X'+str(I)+','+str(K)+','+str(J)]] = 1.0
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                pb.xupper[ix_['KINF'+str(I)+','+str(S)]] = v_['KFRESH']
        v_['LOuKEFF'] = 0.9
        v_['UPuKEFF'] = 1.5
        v_['TEMP'] = float(v_['T'])
        v_['TEMP'] = -0.015*v_['TEMP']
        v_['LOuKEFF'] = v_['LOuKEFF']+v_['TEMP']
        v_['UPuKEFF'] = v_['UPuKEFF']+v_['TEMP']
        for S in range(int(v_['1']),int(v_['T'])+1):
            pb.xlower[ix_['KEFF'+str(S)]] = v_['LOuKEFF']
            pb.xupper[ix_['KEFF'+str(S)]] = v_['UPuKEFF']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        v_['R14'] = 14.0
        v_['TEMP'] = v_['KFRESH']/v_['R14']
        for S in range(int(v_['1']),int(v_['T'])+1):
            if('KEFF'+str(S) in ix_):
                pb.x0[ix_['KEFF'+str(S)]] = float(v_['KEFFuINI'])
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KEFF'+str(S)]),float(v_['KEFFuINI'])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                if('KINF'+str(I)+','+str(S) in ix_):
                    pb.x0[ix_['KINF'+str(I)+','+str(S)]] = float(v_['KFRESH'])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KINF'+str(I)+','+str(S)]),float(v_['KFRESH'])))
                if('PHI'+str(I)+','+str(S) in ix_):
                    pb.x0[ix_['PHI'+str(I)+','+str(S)]] = float(v_['TEMP'])
                else:
                    pb.y0  = (
                          arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['PHI'+str(I)+','+str(S)]),float(v_['TEMP'])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for K in range(int(v_['1']),int(v_['L'])+1):
                for J in range(int(v_['1']),int(v_['M'])+1):
                    if('X'+str(I)+','+str(K)+','+str(J) in ix_):
                        pb.x0[ix_['X'+str(I)+','+str(K)+','+str(J)]] = float(0.5)
                    else:
                        pb.y0  = (
                              arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(K)+','+str(J)]),float(0.5)))
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
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        v_['K'] = 2
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Au'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype, ie, iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 3
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Bu'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype, ie, iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        v_['K'] = 4
        v_['K1'] = -1+v_['K']
        for J in range(int(v_['1']),int(v_['M'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ename = 'Cu'+str(I)+','+str(II)+','+str(J)
                    [ie,ie_,_] = s2mpj_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'en3PROD')
                    ielftype = arrset(ielftype, ie, iet_["en3PROD"])
                    vname = 'X'+str(I)+','+str(int(v_['K']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'X'+str(II)+','+str(int(v_['K1']))+','+str(J)
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                    vname = 'KINF'+str(II)+','+str(int(v_['T']))
                    [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='V3')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ename = 'KTP'+str(I)+','+str(S)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype, ie, iet_["en2PROD"])
                vname = 'KEFF'+str(S)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PHI'+str(I)+','+str(S)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ename = 'P'+str(I)+','+str(S)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'en2PROD')
                ielftype = arrset(ielftype, ie, iet_["en2PROD"])
                vname = 'KINF'+str(I)+','+str(S)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PHI'+str(I)+','+str(S)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['M'])+1):
                for II in range(int(v_['1']),int(v_['N'])+1):
                    ig = ig_['PLAC'+str(I)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['Au'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['V'+str(II)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['Bu'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['V'+str(II)]))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['Cu'+str(I)+','+str(II)+','+str(J)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['V'+str(II)]))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['KERN'+str(I)+','+str(S)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['KTP'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(-1.0))
                v_['TEMP'] = v_['NROW'+str(I)]
                v_['NuROW'] = int(np.fix(v_['TEMP']))
                for II in range(int(v_['1']),int(v_['NuROW'])+1):
                    v_['TEMP'] = v_['ROWu'+str(I)+','+str(II)]
                    v_['III'] = int(np.fix(v_['TEMP']))
                    ig = ig_['KERN'+str(I)+','+str(S)]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['P'+str(int(v_['III']))+','+str(S)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw  = (
                          loaset(pbm.grelw,ig,posel,float(v_['G'+str(I)+','+str(int(v_['III']))])))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T-1'])+1):
                ig = ig_['KINFF'+str(I)+','+str(S)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-ACC']))
        for S in range(int(v_['1']),int(v_['T'])+1):
            for I in range(int(v_['1']),int(v_['N'])+1):
                ig = ig_['CPOW'+str(S)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['V'+str(I)]))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for S in range(int(v_['1']),int(v_['T'])+1):
                ig = ig_['PEAK'+str(I)+','+str(S)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['P'+str(I)+','+str(S)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR2-MN-342-284"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PROD(pbm,nargout,*args):

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
    def en3PROD(pbm,nargout,*args):

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

