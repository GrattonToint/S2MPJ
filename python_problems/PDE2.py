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
#    classification = "LLR2-AN-V-V"
# 
#    the x-y discretization 
# 
#           Alternative values for the SIF file parameters:
# IE N                   3              $-PARAMETER
# IE N                   299            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'PDE2'

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
        if nargin<1:
            v_['N'] = int(29);  #  SIF file default value
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
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['N1'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [iv,ix_,_] = s2mpj_ii('T'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'T'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['0']),int(v_['N1'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['T'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            v_['I+'] = 1+I
            v_['I-'] = -1+I
            for J in range(int(v_['1']),int(v_['N'])+1):
                v_['J+'] = 1+J
                v_['J-'] = -1+J
                [ig,ig_,_] = s2mpj_ii('P'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'P'+str(I)+','+str(J))
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(4.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(I)+','+str(int(v_['J+']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(I)+','+str(int(v_['J-']))]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I+']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(int(v_['I-']))+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'>=')
                cnames = arrset(cnames,ig,'A'+str(I)+','+str(J))
                iv = ix_['T'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['H'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'B'+str(I)+','+str(J))
                iv = ix_['T'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
                iv = ix_['X'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['H'])+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['0'])))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['0'])))
            iv = ix_['X'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['0'])))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['0'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['0'])))
            iv = ix_['X'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['N1'])))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('C'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'C'+str(I)+','+str(int(v_['N1'])))
            iv = ix_['X'+str(I)+','+str(int(v_['N1']))]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['N1'])))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('D'+str(I)+','+str(int(v_['N1'])),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'D'+str(I)+','+str(int(v_['N1'])))
            iv = ix_['X'+str(I)+','+str(int(v_['N1']))]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['0']))+','+str(I))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['0']))+','+str(I))
            iv = ix_['X'+str(int(v_['0']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['0']))+','+str(I))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['0']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['0']))+','+str(I))
            iv = ix_['X'+str(int(v_['0']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['N1']))+','+str(I))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('E'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'E'+str(int(v_['N1']))+','+str(I))
            iv = ix_['X'+str(int(v_['N1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['N1']))+','+str(I))
            iv = ix_['T'+str(I)+','+str(int(v_['0']))]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('F'+str(int(v_['N1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'F'+str(int(v_['N1']))+','+str(I))
            iv = ix_['X'+str(int(v_['N1']))+','+str(I)]
            pbm.A[ig,iv] = float(v_['SQRTAH'])+pbm.A[ig,iv]
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
        for I in range(int(v_['1']),int(v_['N'])+1):
            for J in range(int(v_['1']),int(v_['N'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['P'+str(I)+','+str(J)],float(v_['GH2']))
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
                pbm.gconst = arrset(pbm.gconst,ig_['A'+str(I)+','+str(J)],float(v_['YD']))
                pbm.gconst = arrset(pbm.gconst,ig_['B'+str(I)+','+str(J)],float(v_['-YD']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),3.5)
        for I in range(int(v_['1']),int(v_['N'])+1):
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['0']))]] = 10.0
            pb.xupper[ix_['X'+str(I)+','+str(int(v_['N1']))]] = 10.0
            pb.xupper[ix_['X'+str(int(v_['0']))+','+str(I)]] = 10.0
            pb.xupper[ix_['X'+str(int(v_['N1']))+','+str(I)]] = 10.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "LLR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

