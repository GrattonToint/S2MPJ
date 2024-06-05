from s2mpjlib import *
class  FERRISDC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FERRISDC
#    *********
# 
#    A QP suggested by Michael Ferris
#    classification = ""
#    SIF input: Nick Gould, November 2001.
# 
#    classification = "QLR2-AN-V-V"
# 
#           Alternative values for the SIF file parameters:
# IE n                   4              $-PARAMETER
# IE n                   100            $-PARAMETER
# IE n                   200            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FERRISDC'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'FERRISDC'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['n'] = int(4);  #  SIF file default value
        else:
            v_['n'] = int(args[0])
# IE n                   300            $-PARAMETER
# IE k                   3              $-PARAMETER
# IE k                   10             $-PARAMETER
        if nargin<2:
            v_['k'] = int(3);  #  SIF file default value
        else:
            v_['k'] = int(args[1])
# IE k                   20             $-PARAMETER
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['12'] = 12.0
        v_['24'] = 24.0
        v_['240'] = 240.0
        v_['k-1'] = -1+v_['k']
        v_['k'] = float(v_['k'])
        v_['k-1'] = float(v_['k-1'])
        v_['-1/k-1'] = -1.0/v_['k-1']
        v_['-1/k'] = -1.0/v_['k']
        v_['n'] = float(v_['n'])
        v_['1'] = 1.0
        v_['2'] = 2.0
        v_['1/12'] = v_['1']/v_['12']
        v_['1/24'] = v_['1']/v_['24']
        v_['1/240'] = v_['1']/v_['240']
        v_['7/240'] = 7.0*v_['1/240']
        v_['2**2'] = v_['2']*v_['2']
        v_['2**4'] = v_['2**2']*v_['2**2']
        v_['2**8'] = v_['2**4']*v_['2**4']
        v_['2**10'] = v_['2**8']*v_['2**2']
        v_['2**16'] = v_['2**8']*v_['2**8']
        v_['2**26'] = v_['2**16']*v_['2**10']
        v_['2**-26'] = v_['1']/v_['2**26']
        v_['nlambda'] = v_['n']*v_['2**-26']
        v_['-1/k-1*nl'] = v_['nlambda']*v_['-1/k-1']
        v_['ix'] = 1
        v_['ax'] = 16807
        v_['b15'] = 32768
        v_['b16'] = 65536
        v_['pp'] = 2147483647
        v_['pp'] = float(v_['pp'])
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['xhi'] = int(np.fix(v_['ix']/v_['b16']))
            v_['xalo'] = v_['xhi']*v_['b16']
            v_['xalo'] = v_['ix']-v_['xalo']
            v_['xalo'] = v_['xalo']*v_['ax']
            v_['leftlo'] = int(np.fix(v_['xalo']/v_['b16']))
            v_['fhi'] = v_['xhi']*v_['ax']
            v_['fhi'] = v_['fhi']+v_['leftlo']
            v_['kk'] = int(np.fix(v_['fhi']/v_['b15']))
            v_['dum'] = v_['leftlo']*v_['b16']
            v_['dum'] = v_['xalo']-v_['dum']
            v_['ix'] = v_['dum']-v_['pp']
            v_['dum'] = v_['kk']*v_['b15']
            v_['dum'] = v_['fhi']-v_['dum']
            v_['dum'] = v_['dum']*v_['b16']
            v_['ix'] = v_['ix']+v_['dum']
            v_['ix'] = v_['ix']+v_['kk']
            v_['a'] = float(v_['ix'])
            v_['a'] = -1.0*v_['a']
            v_['b'] = 0.0
            v_['absa'] = np.absolute(v_['a'])
            v_['absb'] = np.absolute(v_['b'])
            v_['absa+b'] = v_['absa']+v_['absb']
            v_['absa+b+2'] = 2.0+v_['absa+b']
            v_['a'] = v_['a']+v_['absa+b+2']
            v_['b'] = v_['b']+v_['absa+b+2']
            v_['a/b'] = v_['a']/v_['b']
            v_['b/a'] = v_['b']/v_['a']
            v_['a/b'] = int(np.fix(v_['a/b']))
            v_['b/a'] = int(np.fix(v_['b/a']))
            v_['a/b'] = float(v_['a/b'])
            v_['b/a'] = float(v_['b/a'])
            v_['sum'] = v_['a/b']+v_['b/a']
            v_['a'] = v_['a']*v_['a/b']
            v_['b'] = v_['b']*v_['b/a']
            v_['maxa,b'] = v_['a']+v_['b']
            v_['maxa,b'] = v_['maxa,b']/v_['sum']
            v_['c'] = v_['absa+b+2']-v_['maxa,b']
            v_['a'] = v_['absa+b+2']-v_['a']
            v_['absc'] = np.absolute(v_['c'])
            v_['absc+1'] = 1.0+v_['absc']
            v_['absc+2'] = 2.0+v_['absc']
            v_['f'] = v_['absc+2']/v_['absc+1']
            v_['f'] = int(np.fix(v_['f']))
            v_['g'] = 2-v_['f']
            for l in range(int(v_['1']),int(v_['g'])+1):
                v_['ix'] = v_['ix']+v_['pp']
            v_['randp'] = float(v_['ix'])
            v_['X'+str(j)] = v_['randp']/v_['pp']
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['xhi'] = int(np.fix(v_['ix']/v_['b16']))
            v_['xalo'] = v_['xhi']*v_['b16']
            v_['xalo'] = v_['ix']-v_['xalo']
            v_['xalo'] = v_['xalo']*v_['ax']
            v_['leftlo'] = int(np.fix(v_['xalo']/v_['b16']))
            v_['fhi'] = v_['xhi']*v_['ax']
            v_['fhi'] = v_['fhi']+v_['leftlo']
            v_['kk'] = int(np.fix(v_['fhi']/v_['b15']))
            v_['dum'] = v_['leftlo']*v_['b16']
            v_['dum'] = v_['xalo']-v_['dum']
            v_['ix'] = v_['dum']-v_['pp']
            v_['dum'] = v_['kk']*v_['b15']
            v_['dum'] = v_['fhi']-v_['dum']
            v_['dum'] = v_['dum']*v_['b16']
            v_['ix'] = v_['ix']+v_['dum']
            v_['ix'] = v_['ix']+v_['kk']
            v_['a'] = float(v_['ix'])
            v_['a'] = -1.0*v_['a']
            v_['b'] = 0.0
            v_['absa'] = np.absolute(v_['a'])
            v_['absb'] = np.absolute(v_['b'])
            v_['absa+b'] = v_['absa']+v_['absb']
            v_['absa+b+2'] = 2.0+v_['absa+b']
            v_['a'] = v_['a']+v_['absa+b+2']
            v_['b'] = v_['b']+v_['absa+b+2']
            v_['a/b'] = v_['a']/v_['b']
            v_['b/a'] = v_['b']/v_['a']
            v_['a/b'] = int(np.fix(v_['a/b']))
            v_['b/a'] = int(np.fix(v_['b/a']))
            v_['a/b'] = float(v_['a/b'])
            v_['b/a'] = float(v_['b/a'])
            v_['sum'] = v_['a/b']+v_['b/a']
            v_['a'] = v_['a']*v_['a/b']
            v_['b'] = v_['b']*v_['b/a']
            v_['maxa,b'] = v_['a']+v_['b']
            v_['maxa,b'] = v_['maxa,b']/v_['sum']
            v_['c'] = v_['absa+b+2']-v_['maxa,b']
            v_['a'] = v_['absa+b+2']-v_['a']
            v_['absc'] = np.absolute(v_['c'])
            v_['absc+1'] = 1.0+v_['absc']
            v_['absc+2'] = 2.0+v_['absc']
            v_['f'] = v_['absc+2']/v_['absc+1']
            v_['f'] = int(np.fix(v_['f']))
            v_['g'] = 2-v_['f']
            for l in range(int(v_['1']),int(v_['g'])+1):
                v_['ix'] = v_['ix']+v_['pp']
            v_['randp'] = float(v_['ix'])
            v_['R'+str(j)] = v_['randp']/v_['pp']
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['arg'] = -3.0*v_['X'+str(j)]
            v_['arg'] = np.exp(v_['arg'])
            v_['P'+str(j)+','+str(int(v_['1']))] = 0.97*v_['arg']
            v_['arg'] = -1.2+v_['X'+str(j)]
            v_['arg'] = v_['arg']*v_['arg']
            v_['arg'] = -2.5*v_['arg']
            v_['P'+str(j)+','+str(int(v_['3']))] = np.exp(v_['arg'])
            v_['arg'] = v_['1']-v_['P'+str(j)+','+str(int(v_['1']))]
            v_['P'+str(j)+','+str(int(v_['2']))] = (v_['arg']-v_['P'+str(j)+','+
                 str(int(v_['3']))])
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['a'] = v_['P'+str(j)+','+str(int(v_['1']))]-v_['R'+str(j)]
            v_['a'] = -1.0*v_['a']
            v_['b'] = 0.0
            v_['absa'] = np.absolute(v_['a'])
            v_['absb'] = np.absolute(v_['b'])
            v_['absa+b'] = v_['absa']+v_['absb']
            v_['absa+b+2'] = 2.0+v_['absa+b']
            v_['a'] = v_['a']+v_['absa+b+2']
            v_['b'] = v_['b']+v_['absa+b+2']
            v_['a/b'] = v_['a']/v_['b']
            v_['b/a'] = v_['b']/v_['a']
            v_['a/b'] = int(np.fix(v_['a/b']))
            v_['b/a'] = int(np.fix(v_['b/a']))
            v_['a/b'] = float(v_['a/b'])
            v_['b/a'] = float(v_['b/a'])
            v_['sum'] = v_['a/b']+v_['b/a']
            v_['a'] = v_['a']*v_['a/b']
            v_['b'] = v_['b']*v_['b/a']
            v_['maxa,b'] = v_['a']+v_['b']
            v_['maxa,b'] = v_['maxa,b']/v_['sum']
            v_['c'] = v_['absa+b+2']-v_['maxa,b']
            v_['a'] = v_['absa+b+2']-v_['a']
            v_['absc'] = np.absolute(v_['c'])
            v_['absc+1'] = 1.0+v_['absc']
            v_['absc+2'] = 2.0+v_['absc']
            v_['f'] = v_['absc+2']/v_['absc+1']
            v_['f'] = int(np.fix(v_['f']))
            v_['g'] = 2-v_['f']
            for l1 in range(int(v_['g']),int(v_['0'])+1):
                v_['y'+str(j)] = 1.0
            for l1 in range(int(v_['1']),int(v_['g'])+1):
                v_['a'] = v_['1']-v_['P'+str(j)+','+str(int(v_['3']))]
                v_['a'] = v_['a']-v_['R'+str(j)]
                v_['a'] = -1.0*v_['a']
                v_['b'] = 0.0
                v_['absa'] = np.absolute(v_['a'])
                v_['absb'] = np.absolute(v_['b'])
                v_['absa+b'] = v_['absa']+v_['absb']
                v_['absa+b+2'] = 2.0+v_['absa+b']
                v_['a'] = v_['a']+v_['absa+b+2']
                v_['b'] = v_['b']+v_['absa+b+2']
                v_['a/b'] = v_['a']/v_['b']
                v_['b/a'] = v_['b']/v_['a']
                v_['a/b'] = int(np.fix(v_['a/b']))
                v_['b/a'] = int(np.fix(v_['b/a']))
                v_['a/b'] = float(v_['a/b'])
                v_['b/a'] = float(v_['b/a'])
                v_['sum'] = v_['a/b']+v_['b/a']
                v_['a'] = v_['a']*v_['a/b']
                v_['b'] = v_['b']*v_['b/a']
                v_['maxa,b'] = v_['a']+v_['b']
                v_['maxa,b'] = v_['maxa,b']/v_['sum']
                v_['c'] = v_['absa+b+2']-v_['maxa,b']
                v_['a'] = v_['absa+b+2']-v_['a']
                v_['absc'] = np.absolute(v_['c'])
                v_['absc+1'] = 1.0+v_['absc']
                v_['absc+2'] = 2.0+v_['absc']
                v_['f'] = v_['absc+2']/v_['absc+1']
                v_['f'] = int(np.fix(v_['f']))
                v_['g'] = 2-v_['f']
                for l2 in range(int(v_['g']),int(v_['0'])+1):
                    v_['y'+str(j)] = 2.0
                for l2 in range(int(v_['1']),int(v_['g'])+1):
                    v_['y'+str(j)] = 3.0
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['yj'] = v_['y'+str(j)]
            v_['yj'] = int(np.fix(v_['yj']))
            for i in range(int(v_['1']),int(v_['k'])+1):
                v_['c'] = v_['yj']-i
                v_['c'] = float(v_['c'])
                v_['absc'] = np.absolute(v_['c'])
                v_['absc+1'] = 1.0+v_['absc']
                v_['absc+2'] = 2.0+v_['absc']
                v_['f'] = v_['absc+2']/v_['absc+1']
                v_['f'] = int(np.fix(v_['f']))
                v_['g'] = 2-v_['f']
                for l in range(int(v_['g']),int(v_['0'])+1):
                    v_['Y'+str(i)+','+str(j)] = v_['nlambda']
                for l in range(int(v_['1']),int(v_['g'])+1):
                    v_['Y'+str(i)+','+str(j)] = v_['-1/k-1*nl']
        for i in range(int(v_['1']),int(v_['n'])+1):
            v_['di'] = -0.5+v_['X'+str(i)]
            v_['di2'] = v_['di']*v_['di']
            v_['di2'] = v_['di2']-v_['1/12']
            for j in range(int(i),int(v_['n'])+1):
                v_['Xi-Xj'] = v_['X'+str(i)]-v_['X'+str(j)]
                v_['bij'] = np.absolute(v_['Xi-Xj'])
                v_['dj'] = -0.5+v_['X'+str(j)]
                v_['dj2'] = v_['dj']*v_['dj']
                v_['dj2'] = v_['dj2']-v_['1/12']
                v_['c'] = -0.5+v_['bij']
                v_['c2'] = v_['c']*v_['c']
                v_['c4'] = v_['c2']*v_['c2']
                v_['c2'] = -0.5*v_['c2']
                v_['arg'] = v_['7/240']+v_['c2']
                v_['arg'] = v_['arg']+v_['c4']
                v_['arg'] = v_['arg']*v_['1/24']
                v_['dij'] = v_['di']*v_['dj']
                v_['dij2'] = v_['di2']*v_['dj2']
                v_['dij2'] = 0.25*v_['dij2']
                v_['arg'] = v_['dij2']-v_['arg']
                v_['K'+str(i)+','+str(j)] = v_['dij']+v_['arg']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for i in range(int(v_['1']),int(v_['k'])+1):
            for j in range(int(v_['1']),int(v_['n'])+1):
                [iv,ix_,_] = s2mpj_ii('A'+str(i)+','+str(j),ix_)
                pb.xnames=arrset(pb.xnames,iv,'A'+str(i)+','+str(j))
        for i in range(int(v_['1']),int(v_['n'])+1):
            [iv,ix_,_] = s2mpj_ii('W'+str(i),ix_)
            pb.xnames=arrset(pb.xnames,iv,'W'+str(i))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for j in range(int(v_['1']),int(v_['n'])+1):
            for i in range(int(v_['1']),int(v_['k'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
                gtype = arrset(gtype,ig,'<>')
                iv = ix_['A'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(v_['Y'+str(i)+','+str(j)])+pbm.A[ig,iv]
        for i in range(int(v_['1']),int(v_['k'])+1):
            for j in range(int(v_['1']),int(v_['n'])+1):
                [ig,ig_,_] = s2mpj_ii('C'+str(i),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'C'+str(i))
                iv = ix_['A'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
                iv = ix_['W'+str(j)]
                pbm.A[ig,iv] = float(v_['-1/k'])+pbm.A[ig,iv]
        for j in range(int(v_['1']),int(v_['n'])+1):
            [ig,ig_,_] = s2mpj_ii('A'+str(j),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'A'+str(j))
            iv = ix_['W'+str(j)]
            pbm.A[ig,iv] = float(-1.0)+pbm.A[ig,iv]
            for i in range(int(v_['1']),int(v_['k'])+1):
                [ig,ig_,_] = s2mpj_ii('A'+str(j),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'A'+str(j))
                iv = ix_['A'+str(i)+','+str(j)]
                pbm.A[ig,iv] = float(1.0)+pbm.A[ig,iv]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),0.0)
        pb.xupper = np.full((pb.n,1),1.0)
        for i in range(int(v_['1']),int(v_['n'])+1):
            pb.xlower[ix_['W'+str(i)]] = -float('Inf')
            pb.xupper[ix_['W'+str(i)]] = +float('Inf')
        for j in range(int(v_['1']),int(v_['n'])+1):
            v_['yj'] = v_['y'+str(j)]
            v_['yj'] = int(np.fix(v_['yj']))
            for i in range(int(v_['1']),int(v_['k'])+1):
                v_['c'] = v_['yj']-i
                v_['c'] = float(v_['c'])
                v_['absc'] = np.absolute(v_['c'])
                v_['absc+1'] = 1.0+v_['absc']
                v_['absc+2'] = 2.0+v_['absc']
                v_['f'] = v_['absc+2']/v_['absc+1']
                v_['f'] = int(np.fix(v_['f']))
                v_['g'] = 2-v_['f']
                for l in range(int(v_['g']),int(v_['0'])+1):
                    pb.xlower[ix_['A'+str(i)+','+str(j)]] = 0.0
                    pb.xupper[ix_['A'+str(i)+','+str(j)]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = lil_matrix((pb.n, pb.n))
        for i in range(int(v_['1']),int(v_['k'])+1):
            for l in range(int(v_['1']),int(v_['n'])+1):
                for j in range(int(v_['1']),int(l)+1):
                    ix1 = ix_['A'+str(i)+','+str(j)]
                    ix2 = ix_['A'+str(i)+','+str(l)]
                    pbm.H[ix1,ix2] = float(v_['K'+str(j)+','+str(l)])+pbm.H[ix1,ix2]
                    pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        for l in range(int(v_['1']),int(v_['n'])+1):
            for j in range(int(v_['1']),int(l)+1):
                v_['coef'] = v_['-1/k']*v_['K'+str(j)+','+str(l)]
                ix1 = ix_['W'+str(j)]
                ix2 = ix_['W'+str(l)]
                pbm.H[ix1,ix2] = float(v_['coef'])+pbm.H[ix1,ix2]
                pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -1.131846D+2   $ nlambda = 1.5625
# XL SOLUTION            -8.032841E-5   $ nlambda = 1.4901E-06
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
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.lincons   = np.arange(len(pbm.congrps))
        pb.pbclass = "QLR2-AN-V-V"
        pb.x0          = np.zeros((pb.n,1))
        pbm.H = pbm.H.tocsr()
        self.pb = pb; self.pbm = pbm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

