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
#    classification = "C-CQOR2-RN-90-78"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'LAKES'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
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
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['5'])+1):
                [iv,ix_,_] = s2mpj_ii('T'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'T'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('O'+str(I)+','+str(J),ix_)
                self.xnames=arrset(self.xnames,iv,'O'+str(I)+','+str(J))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['2'])),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(I)+','+str(int(v_['2'])))
            [iv,ix_,_] = s2mpj_ii('W'+str(I)+','+str(int(v_['2'])),ix_)
            self.xnames=arrset(self.xnames,iv,'W'+str(I)+','+str(int(v_['2'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['3'])),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(I)+','+str(int(v_['3'])))
            [iv,ix_,_] = s2mpj_ii('W'+str(I)+','+str(int(v_['3'])),ix_)
            self.xnames=arrset(self.xnames,iv,'W'+str(I)+','+str(int(v_['3'])))
            [iv,ix_,_] = s2mpj_ii('V'+str(I)+','+str(int(v_['4'])),ix_)
            self.xnames=arrset(self.xnames,iv,'V'+str(I)+','+str(int(v_['4'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for i in range(int(v_['1']),int(v_['N'])+1):
            v_['n+i'] = v_['N']+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                [ig,ig_,_] = s2mpj_ii('R'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(i)+','+str(j)]])
                valA = np.append(valA,float(1.0))
                [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['n+i']))+','+str(j),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['O'+str(i)+','+str(j)]])
                valA = np.append(valA,float(1.0))
        for i in range(int(v_['1']),int(v_['N-1'])+1):
            v_['k+1'] = 1+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                [ig,ig_,_] = s2mpj_ii('G'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(i)+','+str(j))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(int(v_['k+1']))+','+str(j)]])
                valA = np.append(valA,float(1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['T'+str(i)+','+str(j)]])
                valA = np.append(valA,float(-1.0))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['O'+str(i)+','+str(j)]])
                valA = np.append(valA,float(1.0))
            for j in range(int(v_['2']),int(v_['5'])+1):
                v_['j-1'] = -1+j
                [ig,ig_,_] = s2mpj_ii('G'+str(i)+','+str(j),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'G'+str(i)+','+str(j))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['O'+str(i)+','+str(int(v_['j-1']))]])
                valA = np.append(valA,float(-1.0))
        for j in range(int(v_['1']),int(v_['5'])+1):
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(int(v_['1']))+','+str(j)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(int(v_['N']))+','+str(j)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['O'+str(int(v_['N']))+','+str(j)]])
            valA = np.append(valA,float(1.0))
        for j in range(int(v_['2']),int(v_['5'])+1):
            v_['j-1'] = -1+j
            [ig,ig_,_] = s2mpj_ii('G'+str(int(v_['N']))+','+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'G'+str(int(v_['N']))+','+str(j))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['O'+str(int(v_['N']))+','+str(int(v_['j-1']))]])
            valA = np.append(valA,float(-1.0))
        for i in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['O'+str(i)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['O'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('A'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['O'+str(i)+','+str(int(v_['4']))]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(i)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(-1.0))
            v_['C'] = 961.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(v_['C']))
            v_['C'] = 9.2
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(v_['C']))
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(i)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(-1.0))
            v_['C'] = 480.8
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['2']))]])
            valA = np.append(valA,float(v_['C']))
            v_['C'] = -4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['2'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['2'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(v_['C']))
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(-1.0))
            v_['C'] = 4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(v_['C']))
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['W'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(-1.0))
            v_['C'] = 4.6
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['3']))]])
            valA = np.append(valA,float(v_['C']))
            v_['C'] = -105.15
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('W'+str(i)+','+str(int(v_['3'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'W'+str(i)+','+str(int(v_['3'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['4']))]])
            valA = np.append(valA,float(v_['C']))
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['4'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['V'+str(i)+','+str(int(v_['4']))]])
            valA = np.append(valA,float(-1.0))
            v_['C'] = 105.15
            v_['C'] = 1.0/v_['C']
            [ig,ig_,_] = s2mpj_ii('V'+str(i)+','+str(int(v_['4'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(i)+','+str(int(v_['4'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['T'+str(i)+','+str(int(v_['4']))]])
            valA = np.append(valA,float(v_['C']))
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
              arrset(self.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['1']))],float(v_['S1s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['2']))],float(v_['S1s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['3']))],float(v_['S1s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['4']))],float(v_['S1s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['1']))+','+str(int(v_['5']))],float(v_['S1s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['1']))],float(v_['S2s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['2']))],float(v_['S2s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['3']))],float(v_['S2s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['4']))],float(v_['S2s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['2']))+','+str(int(v_['5']))],float(v_['S2s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['1']))],float(v_['S3s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['2']))],float(v_['S3s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['3']))],float(v_['S3s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['4']))],float(v_['S3s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['3']))+','+str(int(v_['5']))],float(v_['S3s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['1']))],float(v_['S4s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['2']))],float(v_['S4s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['3']))],float(v_['S4s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['4']))],float(v_['S4s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['4']))+','+str(int(v_['5']))],float(v_['S4s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['1']))],float(v_['S5s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['2']))],float(v_['S5s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['3']))],float(v_['S5s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['4']))],float(v_['S5s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['5']))+','+str(int(v_['5']))],float(v_['S5s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['1']))],float(v_['S6s1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['2']))],float(v_['S6s2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['3']))],float(v_['S6s3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['4']))],float(v_['S6s4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['6']))+','+str(int(v_['5']))],float(v_['S6s5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['1']))],float(v_['O1o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['2']))],float(v_['O1o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['3']))],float(v_['O1o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['4']))],float(v_['O1o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['7']))+','+str(int(v_['5']))],float(v_['O1o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['1']))],float(v_['O2o1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['2']))],float(v_['O2o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['3']))],float(v_['O2o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['4']))],float(v_['O2o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['8']))+','+str(int(v_['5']))],float(v_['O2o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['1']))],float(v_['O3o1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['2']))],float(v_['O3o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['3']))],float(v_['O3o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['4']))],float(v_['O3o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['9']))+','+str(int(v_['5']))],float(v_['O3o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['1']))],float(v_['O4o1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['2']))],float(v_['O4o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['3']))],float(v_['O4o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['4']))],float(v_['O4o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['10']))+','+str(int(v_['5']))],float(v_['O4o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['1']))],float(v_['O5o1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['2']))],float(v_['O5o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['3']))],float(v_['O5o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['4']))],float(v_['O5o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['11']))+','+str(int(v_['5']))],float(v_['O5o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['1']))],float(v_['O6o1'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['2']))],float(v_['O6o2'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['3']))],float(v_['O6o3'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['4']))],float(v_['O6o4'])))
        self.gconst  = (
              arrset(self.gconst,ig_['R'+str(int(v_['12']))+','+str(int(v_['5']))],float(v_['O6o5'])))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['1']))],float(-22.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['2']))],float(-1.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['3']))],float(3.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['4']))],float(-27.2)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['1']))+','+str(int(v_['5']))],float(51.5)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['1']))],float(44.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['2']))],float(162.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['3']))],float(8.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['4']))],float(12.5)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['2']))+','+str(int(v_['5']))],float(53.5)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['1']))],float(-11.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['2']))],float(60.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['3']))],float(10.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['4']))],float(18.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['3']))+','+str(int(v_['5']))],float(39.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['1']))],float(124.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['2']))],float(246.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['3']))],float(6.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['4']))],float(9.7)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['4']))+','+str(int(v_['5']))],float(17.2)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['1']))],float(127.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['2']))],float(175.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['3']))],float(3.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['4']))],float(10.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['5']))+','+str(int(v_['5']))],float(30.2)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['1']))],float(78.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['2']))],float(156.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['3']))],float(3.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['4']))],float(14.0)))
        self.gconst  = (
              arrset(self.gconst,ig_['G'+str(int(v_['6']))+','+str(int(v_['5']))],float(23.2)))
        for i in range(int(v_['1']),int(v_['N'])+1):
            self.gconst  = (
                  arrset(self.gconst,ig_['V'+str(i)+','+str(int(v_['2']))],float(543.4)))
            self.gconst  = (
                  arrset(self.gconst,ig_['W'+str(i)+','+str(int(v_['2']))],float(0.0)))
            self.gconst  = (
                  arrset(self.gconst,ig_['V'+str(i)+','+str(int(v_['3']))],float(543.4)))
            self.gconst  = (
                  arrset(self.gconst,ig_['W'+str(i)+','+str(int(v_['3']))],float(0.0)))
            self.gconst  = (
                  arrset(self.gconst,ig_['V'+str(i)+','+str(int(v_['4']))],float(550.11)))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        for i in range(int(v_['1']),int(v_['N'])+1):
            self.xlower[ix_['W'+str(i)+','+str(int(v_['2']))]] = 0.0001
            self.xlower[ix_['W'+str(i)+','+str(int(v_['3']))]] = 0.0001
            self.xlower[ix_['V'+str(i)+','+str(int(v_['4']))]] = 0.0001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(1.0))
        self.y0 = np.full((self.m,1),float(1.0))
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
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for i in range(int(v_['1']),int(v_['N'])+1):
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2VAR')
            ielftype = arrset(ielftype,ie,iet_["en2VAR"])
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'W'+str(i)+','+str(int(v_['2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='W')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['1']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(0.0841168))
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2VAR')
            ielftype = arrset(ielftype,ie,iet_["en2VAR"])
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'W'+str(i)+','+str(int(v_['3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='W')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['2']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(0.1280849))
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en1VAR')
            ielftype = arrset(ielftype,ie,iet_["en1VAR"])
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'V'+str(i)+','+str(int(v_['4']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),float(1.0))
            posev = np.where(elftv[ielftype[ie]]=='V')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'B'+str(i)+','+str(int(v_['3']))
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            posep = np.where(elftp[ielftype[ie]]=='P')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(0.2605))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['A'+str(i)+','+str(int(v_['1']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['1']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['A'+str(i)+','+str(int(v_['2']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            ig = ig_['A'+str(i)+','+str(int(v_['3']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['B'+str(i)+','+str(int(v_['3']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for i in range(int(v_['1']),int(v_['N'])+1):
            v_['n+i'] = v_['N']+i
            for j in range(int(v_['1']),int(v_['5'])+1):
                ig = ig_['R'+str(i)+','+str(j)]
                self.grftype = arrset(self.grftype,ig,'gL2')
                ig = ig_['R'+str(int(v_['n+i']))+','+str(j)]
                self.grftype = arrset(self.grftype,ig,'gL2')
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
        self.pbclass   = "C-CQOR2-RN-90-78"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en1VAR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0,0]**2.2
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*2.2*EV_[0,0]**1.2
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = self.elpar[iel_][0]*2.64*EV_[0,0]**0.2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en2VAR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = self.elpar[iel_][0]*EV_[0,0]**2*EV_[1,0]**0.5
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*2.0*EV_[0,0]*EV_[1,0]**0.5
            g_[1] = self.elpar[iel_][0]*0.5*EV_[0,0]**2/EV_[1,0]**0.5
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = self.elpar[iel_][0]*2.0*EV_[1,0]**0.5
                H_[0,1] = self.elpar[iel_][0]*EV_[0,0]/EV_[1,0]**0.5
                H_[1,0] = H_[0,1]
                H_[1,1] = -self.elpar[iel_][0]*0.25*EV_[0,0]**2/EV_[1,0]**1.5
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

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

