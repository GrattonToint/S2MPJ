from s2mpjlib import *
class  LAKES(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A problem of water resource management in Canada, which may be 
#    formulated as
#    Min  SUM   SUM  (T(i,j)- R(i,j))^2 + (O(i,j)-R(N+i,j))^2)
#        i=1,N j=1,5 
#    subject to
#    T(i+1,1)-T(i,1)+O(i,1)        =  G(i,1)
#    T(i+1,2)-T(i,2)-O(i,1)+O(i,2) =  G(i,2)
#    T(i+1,3)-T(i,3)-O(i,2)+O(i,3) =  G(i,3)
#    T(i+1,4)-T(i,4)-O(i,3)+O(i,4) =  G(i,4)
#    T(i+1,5)-T(i,5)-O(i,4)+O(i,5) =  G(i,5) 
#    i=1,N and T(N+1,j) = T(1,j)  for j=1,5
#    O(i,2)-a*((T(i,2)/480.8+T(i,3)/4.6)/2-543.4)^2 * 
#    (T(i,2)/480.8-T(i,3)/4.6)^.5=0
#    O(i,3)-b*((T(i,3)/4.6-543.4)^2*(T(i,3)/4.6-T(i,4)/105.15)^0.5) = 0
#    O(i,4)-c*(T(i,4)/105.15-550.11)^2.2 = 0
#    where T(i,j) and O(i,j) are variables, R(i,j) are given and
#    a=.0841168  b=.1280849 and c=0.2605.
#    Extra variables 
#    
#    v(i,2) = T(i,2) / 961.6 + T(i,3) / 9.2 - 543.4
#    w(i,2) = T(i,2) / 480.8 - T(i,3) / 4.6
#    v(i,3) = T(i,3) / 4.6 - 543.4
#    w(i,3) = T(i,3) / 4.6 - T(i,4) / 105.15
#    v(i,4) = T(i,4) / 105.15 - 550.11
#    are introduced so that the nonlinear constraints may be rewritten as
#    O(i,2)-a*v(i,2)^2 * w(i,2)^0.5 = 0 ; w(i,2) > 0
#    O(i,3)-b*v(i,3)^2 * w(i,3)^0.5 = 0 ; w(i,3) > 0
#    O(i,4)-c*v(i,4)^2.2 = 0 ; v(i,4) > 0
#    Source:
#    S Jafar Sadjadi
#    Dept. of Systems Design Engineering
#    University of Waterloo
#    Ontario, N2L 3G1 Canada
# 
#    SIF input: Nick Gould and Jafar Sadjadi, November 1995
# 
#    classification = "QOR2-RN-90-78"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LAKES'

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
        v_['N'] = 6
        v_['N-1'] = -1+v_['N']
        v_['N+1'] = 1+v_['N']
        v_['NN'] = 2*v_['N']
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
        v_['S1s1'] = 202761.072
        v_['S1s2'] = 277791.816
        v_['S1s3'] = 2636.996
        v_['S1s4'] = 59987.0235
        v_['S1s5'] = 19490.4
        v_['S2s1'] = 202703.646
        v_['S2s2'] = 277849.512
        v_['S2s3'] = 2638.1
        v_['S2s4'] = 59998.59
        v_['S2s5'] = 19555.2
        v_['S3s1'] = 202720.536
        v_['S3s2'] = 277955.288
        v_['S3s3'] = 2639.894
        v_['S3s4'] = 60046.959
        v_['S3s5'] = 19597.6
        v_['S4s1'] = 202808.364
        v_['S4s2'] = 278104.336
        v_['S4s3'] = 2640.906
        v_['S4s4'] = 60074.298
        v_['S4s5'] = 19652.8
        v_['S5s1'] = 202916.46
        v_['S5s2'] = 278224.536
        v_['S5s3'] = 2641.458
        v_['S5s4'] = 60091.122
        v_['S5s5'] = 19708.8
        v_['S6s1'] = 202953.618
        v_['S6s2'] = 278277.424
        v_['S6s3'] = 2641.458
        v_['S6s4'] = 60082.71
        v_['S6s5'] = 19706.4
        v_['O1o1'] = 83.728
        v_['O1o2'] = 174.665
        v_['O1o3'] = 180.539
        v_['O1o4'] = 211.558
        v_['O1o5'] = 232.252
        v_['O2o1'] = 83.789
        v_['O2o2'] = 173.255
        v_['O2o3'] = 179.917
        v_['O2o4'] = 210.585
        v_['O2o5'] = 215.254
        v_['O3o1'] = 82.9160
        v_['O3o2'] = 173.721
        v_['O3o3'] = 182.676
        v_['O3o4'] = 207.838
        v_['O3o5'] = 203.855
        v_['O4o1'] = 80.134
        v_['O4o2'] = 178.654
        v_['O4o3'] = 185.917
        v_['O4o4'] = 206.416
        v_['O4o5'] = 186.308
        v_['O5o1'] = 65.345
        v_['O5o2'] = 188.01
        v_['O5o3'] = 192.568
        v_['O5o4'] = 204.3
        v_['O5o5'] = 201.1
        v_['O6o1'] = 72.005
        v_['O6o2'] = 193.833
        v_['O6o3'] = 196.651
        v_['O6o4'] = 204.25
        v_['O6o5'] = 241.079
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                [iv,ix_,_] = s2mpj_ii('T'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'T'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('O'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'O'+str(I)+','+str(J))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['2'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(int(v_['2'])))
            [iv,ix_,_] = s2mpj_ii('W'+str(I)+','+str(int(v_['2'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I)+','+str(int(v_['2'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['3'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(int(v_['3'])))
            [iv,ix_,_] = s2mpj_ii('W'+str(I)+','+str(int(v_['3'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(I)+','+str(int(v_['3'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['4'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(int(v_['4'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for i in range(int(v_['1']),int(v_['N'])+1):
            v_['n+i'] = v_['N']+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                [ig,ig_,_] = s2mpj_ii('R'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['T'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['n+i']))+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['O'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['N-1'])+1):
            v_['k+1'] = 1+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(i)+','+str(j))
                iv = ix_['T'+str(int(v_['k+1']))+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['T'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['O'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            for j in range(int(v_['2']),int(v_['5'])+1):
                v_['j-1'] = -1+j
                [ig,ig_,_] = s2mpj_ii('G'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(i)+','+str(j))
                iv = ix_['O'+str(i)+','+str(int(v_['j-1']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for j in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            iv = ix_['T'+str(int(v_['1']))+','+str(j)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            iv = ix_['T'+str(int(v_['N']))+','+str(j)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            iv = ix_['O'+str(int(v_['N']))+','+str(j)]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for j in range(int(v_['2']),int(v_['5'])+1):
            v_['j-1'] = -1+j
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            iv = ix_['O'+str(int(v_['N']))+','+str(int(v_['j-1']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['1'])))
            iv = ix_['O'+str(i)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['2'])))
            iv = ix_['O'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['3'])))
            iv = ix_['O'+str(i)+','+str(int(v_['4']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            iv = ix_['V'+str(i)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['C'] = 961.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            iv = ix_['T'+str(i)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            v_['C'] = 9.2
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            iv = ix_['T'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            iv = ix_['W'+str(i)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['C'] = 480.8
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            iv = ix_['T'+str(i)+','+str(int(v_['2']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            v_['C'] = -4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            iv = ix_['T'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['3'])))
            iv = ix_['V'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['C'] = 4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['3'])))
            iv = ix_['T'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            iv = ix_['W'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['C'] = 4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            iv = ix_['T'+str(i)+','+str(int(v_['3']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            v_['C'] = -105.15
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            iv = ix_['T'+str(i)+','+str(int(v_['4']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['4'])))
            iv = ix_['V'+str(i)+','+str(int(v_['4']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            v_['C'] = 105.15
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['4'])))
            iv = ix_['T'+str(i)+','+str(int(v_['4']))]
            pbm.A[ig,iv] = float(v_['C'])+pbm.A[ig,iv]
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
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['1']))],float(v_['S1s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['2']))],float(v_['S1s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['3']))],float(v_['S1s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['4']))],float(v_['S1s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['5']))],float(v_['S1s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['1']))],float(v_['S2s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['2']))],float(v_['S2s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['3']))],float(v_['S2s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['4']))],float(v_['S2s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['5']))],float(v_['S2s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['1']))],float(v_['S3s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['2']))],float(v_['S3s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['3']))],float(v_['S3s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['4']))],float(v_['S3s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['5']))],float(v_['S3s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['1']))],float(v_['S4s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['2']))],float(v_['S4s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['3']))],float(v_['S4s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['4']))],float(v_['S4s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['5']))],float(v_['S4s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['1']))],float(v_['S5s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['2']))],float(v_['S5s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['3']))],float(v_['S5s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['4']))],float(v_['S5s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['5']))],float(v_['S5s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['1']))],float(v_['S6s1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['2']))],float(v_['S6s2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['3']))],float(v_['S6s3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['4']))],float(v_['S6s4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['5']))],float(v_['S6s5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['1']))],float(v_['O1o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['2']))],float(v_['O1o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['3']))],float(v_['O1o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['4']))],float(v_['O1o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['5']))],float(v_['O1o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['1']))],float(v_['O2o1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['2']))],float(v_['O2o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['3']))],float(v_['O2o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['4']))],float(v_['O2o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['5']))],float(v_['O2o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['1']))],float(v_['O3o1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['2']))],float(v_['O3o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['3']))],float(v_['O3o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['4']))],float(v_['O3o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['5']))],float(v_['O3o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['1']))],float(v_['O4o1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['2']))],float(v_['O4o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['3']))],float(v_['O4o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['4']))],float(v_['O4o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['5']))],float(v_['O4o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['1']))],float(v_['O5o1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['2']))],float(v_['O5o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['3']))],float(v_['O5o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['4']))],float(v_['O5o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['5']))],float(v_['O5o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['1']))],float(v_['O6o1'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['2']))],float(v_['O6o2'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['3']))],float(v_['O6o3'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['4']))],float(v_['O6o4'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['5']))],float(v_['O6o5'])))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['1']))],float(-22.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['2']))],float(-1.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['3']))],float(3.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['4']))],float(-27.2)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['5']))],float(51.5)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['1']))],float(44.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['2']))],float(162.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['3']))],float(8.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['4']))],float(12.5)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['5']))],float(53.5)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['1']))],float(-11.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['2']))],float(60.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['3']))],float(10.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['4']))],float(18.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['5']))],float(39.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['1']))],float(124.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['2']))],float(246.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['3']))],float(6.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['4']))],float(9.7)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['5']))],float(17.2)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['1']))],float(127.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['2']))],float(175.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['3']))],float(3.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['4']))],float(10.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['5']))],float(30.2)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['1']))],float(78.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['2']))],float(156.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['3']))],float(3.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['4']))],float(14.0)))
        pbm.gconst  = (
              arrset(pbm.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['5']))],float(23.2)))
        for i in range(int(v_['1']),int(v_['N'])+1):
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['V'+str(i)+','+str(int(v_['2']))],float(543.4)))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['W'+str(i)+','+str(int(v_['2']))],float(0.0)))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['V'+str(i)+','+str(int(v_['3']))],float(543.4)))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['W'+str(i)+','+str(int(v_['3']))],float(0.0)))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['V'+str(i)+','+str(int(v_['4']))],float(550.11)))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        pb.xlower = np.zeros((pb.n,1))
        for i in range(int(v_['1']),int(v_['N'])+1):
            pb.xlower[ix_['W'+str(i)+','+str(int(v_['2']))]] = 0.0001
            pb.xlower[ix_['W'+str(i)+','+str(int(v_['3']))]] = 0.0001
            pb.xlower[ix_['V'+str(i)+','+str(int(v_['4']))]] = 0.0001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(1.0))
        pb.y0 = np.full((pb.m,1),float(1.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en1VAR', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftp = []
        elftp = loaset(elftp,it,0,'P')
        [it,iet_,_] = s2mpj_ii( 'en2VAR', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        elftp = loaset(elftp,it,0,'P')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for i in range(int(v_['1']),int(v_['N'])+1):
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2VAR')
            ielftype = arrset(ielftype, ie, iet_["en2VAR"])
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'W'+str(i)+','+str(int(v_['2']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0841168))
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en2VAR')
            ielftype = arrset(ielftype, ie, iet_["en2VAR"])
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['3']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'W'+str(i)+','+str(int(v_['3']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='W')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.1280849))
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'en1VAR')
            ielftype = arrset(ielftype, ie, iet_["en1VAR"])
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['4']))
            [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,1.0)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = find(elftp[ielftype[ie]],lambda x:x=='P')
            pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.2605))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['A'+str(i)+','+str(int(v_['1']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['A'+str(i)+','+str(int(v_['2']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['2']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['A'+str(i)+','+str(int(v_['3']))]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['3']))])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for i in range(int(v_['1']),int(v_['N'])+1):
            v_['n+i'] = v_['N']+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                ig = ig_['R'+str(i)+','+str(j)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL2')
                ig = ig_['R'+str(int(v_['n+i']))+','+str(j)]
                pbm.grftype = arrset(pbm.grftype,ig,'gL2')
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
        pb.pbclass = "QOR2-RN-90-78"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en1VAR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]**2.2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*2.2*EV_[0]**1.2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = pbm.elpar[iel_][0]*2.64*EV_[0]**0.2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2VAR(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = pbm.elpar[iel_][0]*EV_[0]**2*EV_[1]**0.5
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = pbm.elpar[iel_][0]*2.0*EV_[0]*EV_[1]**0.5
            g_[1] = pbm.elpar[iel_][0]*0.5*EV_[0]**2/EV_[1]**0.5
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = pbm.elpar[iel_][0]*2.0*EV_[1]**0.5
                H_[0,1] = pbm.elpar[iel_][0]*EV_[0]/EV_[1]**0.5
                H_[1,0] = H_[0,1]
                H_[1,1] = -pbm.elpar[iel_][0]*0.25*EV_[0]**2/EV_[1]**1.5
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

