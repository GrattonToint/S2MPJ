from s2mpjlib import *
class  ROTDISC(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Optimal design of a rotating disk of minimal weight, with constraints on
#    stress and profile.
#    The problem arise in the mechanical design of turbine where several disks 
#    are assembled, as in jet engines and steam turbines in power generation
#    systems. The data correspond to the real design problem for the engine of 
#    a small civil jet.
#    The problem has a linear objective function, linear constraints and 
#    quadratic equality constraints.  The solution lies at a vertex of the
#    feasible set.
# 
#    Source:
#    B. Apraxine and E. Loute,
#    "The optimal design of a rotating disk: a test case for nonlinear
#    programming codes",
#    Facultes Universitaires Saint Louis, Brussels, 1993.
#    See also:
#    J. P. Nigoghossian,
#    "Problem: the optimization of jet engine discs",
#    in "Optimisation and Design", M. Avriel, M. J. Rijckaert and D. J. Wilde,
#    eds., Prentice Hall, Englewood Cliffs, 1973.
# 
#    SIF input : E. Loute and Ph. L. Toint, April 1993
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-CLQR2-RN-905-1081"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ROTDISC'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['K'] = 180
        v_['rotspeed'] = 22000.0
        v_['ri'] = 41.0
        v_['ro'] = 131.0
        v_['rim'] = 127.0
        v_['rho'] = 8200.0
        v_['E'] = 18525.0
        v_['epsc'] = 14.0e-6
        v_['nu'] = 0.3
        v_['wro'] = 13.0
        v_['wmin'] = 4.0
        v_['wmax'] = 34.0
        v_['sigmaro'] = 18.5
        v_['sigmari'] = 0.0
        v_['sigmati'] = 72.0
        v_['sigmato'] = 47.5
        v_['sigmatA'] = 66.0
        v_['sigmaru'] = 55.0
        v_['tempo'] = 450.0
        v_['tempn'] = 150.0
        v_['n'] = 4
        v_['ech1'] = 100.0
        v_['ech2'] = 1.0E-13
        v_['ech3'] = 1.0E-9
        v_['1'] = 1
        v_['0'] = 0
        v_['K-1'] = -1+v_['K']
        v_['RK'] = float(v_['K'])
        v_['Dr'] = v_['ro']-v_['ri']
        v_['dr'] = v_['Dr']/v_['RK']
        v_['pi'] = 3.1415926535
        v_['2pi'] = 2.0*v_['pi']
        v_['aux1'] = v_['2pi']*v_['rotspeed']
        v_['60.0'] = 60.0
        v_['omega'] = v_['aux1']/v_['60.0']
        v_['omega2'] = v_['omega']*v_['omega']
        v_['romg2'] = v_['rho']*v_['omega2']
        v_['romg2/2'] = 0.5*v_['romg2']
        v_['aux2'] = v_['romg2/2']*v_['ech2']
        v_['1+nu'] = 1.0+v_['nu']
        v_['3nu'] = 3.0*v_['nu']
        v_['1+3nu'] = 1.0+v_['3nu']
        v_['3+nu'] = 3.0+v_['nu']
        v_['(1+nu)/2'] = 0.5*v_['1+nu']
        v_['(1+3nu)/2'] = 0.5*v_['1+3nu']
        v_['(3+nu)/2'] = 0.5*v_['3+nu']
        v_['pirho'] = v_['pi']*v_['rho']
        v_['dr/2'] = 0.5*v_['dr']
        v_['aux3'] = v_['aux2']*v_['dr']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for k in range(int(v_['0']),int(v_['K'])+1):
            [iv,ix_,_] = s2mpj_ii('w'+str(k),ix_)
            self.xnames=arrset(self.xnames,iv,'w'+str(k))
            [iv,ix_,_] = s2mpj_ii('sigt'+str(k),ix_)
            self.xnames=arrset(self.xnames,iv,'sigt'+str(k))
            [iv,ix_,_] = s2mpj_ii('sigr'+str(k),ix_)
            self.xnames=arrset(self.xnames,iv,'sigr'+str(k))
            [iv,ix_,_] = s2mpj_ii('x'+str(k),ix_)
            self.xnames=arrset(self.xnames,iv,'x'+str(k))
            [iv,ix_,_] = s2mpj_ii('y'+str(k),ix_)
            self.xnames=arrset(self.xnames,iv,'y'+str(k))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        v_['rk'] = v_['ri']
        v_['-dr/2'] = -1.0*v_['dr/2']
        v_['rk2'] = v_['rk']*v_['rk']
        for k in range(int(v_['0']),int(v_['K-1'])+1):
            v_['k+1'] = 1+k
            v_['rk+1'] = v_['rk']+v_['dr']
            v_['coef1'] = v_['aux3']*v_['rk2']
            v_['rk+1sq'] = v_['rk+1']*v_['rk+1']
            v_['coef2'] = v_['aux3']*v_['rk+1sq']
            [ig,ig_,_] = s2mpj_ii('SR'+str(k),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'SR'+str(k))
            self.gscale = arrset(self.gscale,ig,float(v_['ech1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(k)]])
            valA = np.append(valA,float(v_['coef1']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(v_['coef2']))
            v_['tmp1'] = v_['(1+3nu)/2']*v_['rk']
            v_['tmp2'] = v_['(1+nu)/2']*v_['rk+1']
            v_['tmp3'] = v_['tmp1']-v_['tmp2']
            v_['coef3'] = v_['tmp3']/v_['rk']
            [ig,ig_,_] = s2mpj_ii('ST'+str(k),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ST'+str(k))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['sigr'+str(k)]])
            valA = np.append(valA,float(v_['coef3']))
            v_['tmp4'] = v_['(3+nu)/2']*v_['rk']
            v_['tmp5'] = v_['(1+nu)/2']*v_['rk+1']
            v_['tmp6'] = v_['tmp5']-v_['tmp4']
            v_['coef4'] = v_['tmp6']/v_['rk']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['sigt'+str(k)]])
            valA = np.append(valA,float(v_['coef4']))
            v_['tmp7'] = v_['(1+3nu)/2']*v_['rk+1']
            v_['tmp8'] = v_['(1+nu)/2']*v_['rk']
            v_['tmp9'] = v_['tmp8']-v_['tmp7']
            v_['coef5'] = v_['tmp9']/v_['rk+1']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['sigr'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(v_['coef5']))
            v_['tmp10'] = v_['(3+nu)/2']*v_['rk+1']
            v_['tmp11'] = v_['(1+nu)/2']*v_['rk']
            v_['tmp12'] = v_['tmp10']-v_['tmp11']
            v_['coef6'] = v_['tmp12']/v_['rk+1']
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['sigt'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(v_['coef6']))
            [ig,ig_,_] = s2mpj_ii('STAy'+str(k),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'STAy'+str(k))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['y'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['y'+str(k)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('STAx'+str(k),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'STAx'+str(k))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['x'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['x'+str(k)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(v_['-dr/2']))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(k)]])
            valA = np.append(valA,float(v_['-dr/2']))
            v_['rk'] = v_['rk+1']
            v_['rk2'] = v_['rk+1sq']
        v_['rk-1'] = v_['ri']
        v_['rk'] = v_['rk-1']+v_['dr']
        v_['rk-1sq'] = v_['rk-1']*v_['rk-1']
        v_['rk2'] = v_['rk']*v_['rk']
        v_['aux3'] = v_['rk-1sq']-v_['rk2']
        v_['aux4'] = v_['aux3']*v_['pirho']
        v_['coef1'] = v_['aux4']*v_['ech3']
        v_['-coef1'] = -1.0*v_['coef1']
        [ig,ig_,_] = s2mpj_ii('WEIGHT',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['w'+str(int(v_['0']))]])
        valA = np.append(valA,float(v_['-coef1']))
        for k in range(int(v_['1']),int(v_['K-1'])+1):
            v_['k-1'] = -1+k
            v_['rk+1'] = v_['rk']+v_['dr']
            v_['rk+1sq'] = v_['rk+1']*v_['rk+1']
            v_['aux3'] = v_['rk-1sq']-v_['rk+1sq']
            v_['aux4'] = v_['aux3']*v_['pirho']
            v_['coef1'] = v_['aux4']*v_['ech3']
            v_['-coef1'] = -1.0*v_['coef1']
            [ig,ig_,_] = s2mpj_ii('WEIGHT',ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'WEIGHT')
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(k)]])
            valA = np.append(valA,float(v_['-coef1']))
            v_['rk-1sq'] = v_['rk2']
            v_['rk'] = v_['rk+1']
            v_['rk2'] = v_['rk+1sq']
        v_['aux3'] = v_['rk-1sq']-v_['rk2']
        v_['aux4'] = v_['pirho']*v_['aux3']
        v_['coef1'] = v_['aux4']*v_['ech3']
        v_['-coef1'] = -1.0*v_['coef1']
        [ig,ig_,_] = s2mpj_ii('WEIGHT',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['w'+str(int(v_['K']))]])
        valA = np.append(valA,float(v_['-coef1']))
        for k in range(int(v_['0']),int(v_['K-1'])+1):
            v_['k+1'] = 1+k
            [ig,ig_,_] = s2mpj_ii('SLOP'+str(k),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SLOP'+str(k))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(k)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('SLOM'+str(k),ig_)
            gtype = arrset(gtype,ig,'<=')
            cnames = arrset(cnames,ig,'SLOM'+str(k))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(int(v_['k+1']))]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['w'+str(k)]])
            valA = np.append(valA,float(1.0))
        v_['-sigmatA'] = -1.0*v_['sigmatA']
        [ig,ig_,_] = s2mpj_ii('AVsigt',ig_)
        gtype = arrset(gtype,ig,'<=')
        cnames = arrset(cnames,ig,'AVsigt')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['y'+str(int(v_['K']))]])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['x'+str(int(v_['K']))]])
        valA = np.append(valA,float(v_['-sigmatA']))
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
        v_['Eepsc'] = v_['E']*v_['epsc']
        v_['rk'] = v_['ri']
        v_['aux1'] = v_['rk']/v_['ro']
        v_['aux2'] = v_['aux1']*v_['aux1']
        v_['aux2'] = v_['aux2']*v_['aux2']
        v_['tmp1'] = v_['aux2']*v_['tempn']
        v_['tk'] = v_['tmp1']+v_['tempo']
        for k in range(int(v_['0']),int(v_['K-1'])+1):
            v_['rk+1'] = v_['rk']+v_['dr']
            v_['aux1'] = v_['rk+1']/v_['ro']
            v_['aux2'] = v_['aux1']*v_['aux1']
            v_['aux2'] = v_['aux2']*v_['aux2']
            v_['tmp1'] = v_['aux2']*v_['tempn']
            v_['tk+1'] = v_['tmp1']+v_['tempo']
            v_['tmp2'] = v_['tk']-v_['tk+1']
            v_['coef1'] = v_['tmp2']*v_['Eepsc']
            self.gconst = arrset(self.gconst,ig_['ST'+str(k)],float(v_['coef1']))
            v_['rk'] = v_['rk+1']
            v_['tk'] = v_['tk+1']
        v_['4dr'] = 4.0*v_['dr']
        for k in range(int(v_['0']),int(v_['K-1'])+1):
            self.gconst = arrset(self.gconst,ig_['SLOP'+str(k)],float(v_['4dr']))
            self.gconst = arrset(self.gconst,ig_['SLOM'+str(k)],float(v_['4dr']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower[ix_['sigr'+str(int(v_['0']))]] = v_['sigmari']
        self.xupper[ix_['sigr'+str(int(v_['0']))]] = v_['sigmari']
        self.xupper[ix_['sigt'+str(int(v_['0']))]] = v_['sigmati']
        self.xlower[ix_['sigt'+str(int(v_['0']))]] = -100.0
        self.xlower[ix_['x'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['x'+str(int(v_['0']))]] = 0.0
        self.xlower[ix_['y'+str(int(v_['0']))]] = 0.0
        self.xupper[ix_['y'+str(int(v_['0']))]] = 0.0
        for k in range(int(v_['1']),int(v_['K'])+1):
            self.xlower[ix_['x'+str(k)]] = 0.0
        for k in range(int(v_['1']),int(v_['K-1'])+1):
            self.xlower[ix_['sigr'+str(k)]] = 0.0
            self.xupper[ix_['sigr'+str(k)]] = v_['sigmaru']
            self.xlower[ix_['sigt'+str(k)]] = -100.0
            self.xupper[ix_['sigt'+str(k)]] = 100.0
        v_['K-9'] = -9+v_['K']
        for k in range(int(v_['0']),int(v_['K-9'])+1):
            self.xlower[ix_['w'+str(k)]] = v_['wmin']
            self.xupper[ix_['w'+str(k)]] = v_['wmax']
        v_['K-8'] = -8+v_['K']
        for k in range(int(v_['K-8']),int(v_['K'])+1):
            self.xlower[ix_['w'+str(k)]] = v_['wro']
            self.xupper[ix_['w'+str(k)]] = v_['wro']
        self.xupper[ix_['sigt'+str(int(v_['K']))]] = v_['sigmato']
        self.xlower[ix_['sigr'+str(int(v_['K']))]] = v_['sigmaro']
        self.xupper[ix_['sigr'+str(int(v_['K']))]] = v_['sigmaro']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        if('w0' in ix_):
            self.x0[ix_['w0']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w0']),float(30.000000000)))
        if('sigt0' in ix_):
            self.x0[ix_['sigt0']] = float(70.224043123)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt0']),float(70.224043123)))
        if('sigr0' in ix_):
            self.x0[ix_['sigr0']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr0']),float(0.0)))
        if('x0' in ix_):
            self.x0[ix_['x0']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x0']),float(0.0)))
        if('y0' in ix_):
            self.x0[ix_['y0']] = float(0.0)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y0']),float(0.0)))
        if('w1' in ix_):
            self.x0[ix_['w1']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w1']),float(30.000000000)))
        if('sigt1' in ix_):
            self.x0[ix_['sigt1']] = float(69.337178701)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt1']),float(69.337178701)))
        if('sigr1' in ix_):
            self.x0[ix_['sigr1']] = float(.75150203489)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr1']),float(.75150203489)))
        if('x1' in ix_):
            self.x0[ix_['x1']] = float(15.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x1']),float(15.000000000)))
        if('y1' in ix_):
            self.x0[ix_['y1']] = float(1046.7091637)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y1']),float(1046.7091637)))
        if('w2' in ix_):
            self.x0[ix_['w2']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w2']),float(30.000000000)))
        if('sigt2' in ix_):
            self.x0[ix_['sigt2']] = float(68.478656583)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt2']),float(68.478656583)))
        if('sigr2' in ix_):
            self.x0[ix_['sigr2']] = float(1.4725717269)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr2']),float(1.4725717269)))
        if('x2' in ix_):
            self.x0[ix_['x2']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x2']),float(30.000000000)))
        if('y2' in ix_):
            self.x0[ix_['y2']] = float(2080.3279283)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y2']),float(2080.3279283)))
        if('w3' in ix_):
            self.x0[ix_['w3']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w3']),float(30.000000000)))
        if('sigt3' in ix_):
            self.x0[ix_['sigt3']] = float(67.647086003)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt3']),float(67.647086003)))
        if('sigr3' in ix_):
            self.x0[ix_['sigr3']] = float(2.1645828149)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr3']),float(2.1645828149)))
        if('x3' in ix_):
            self.x0[ix_['x3']] = float(45.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x3']),float(45.000000000)))
        if('y3' in ix_):
            self.x0[ix_['y3']] = float(3101.2709977)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y3']),float(3101.2709977)))
        if('w4' in ix_):
            self.x0[ix_['w4']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w4']),float(30.000000000)))
        if('sigt4' in ix_):
            self.x0[ix_['sigt4']] = float(66.841155637)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt4']),float(66.841155637)))
        if('sigr4' in ix_):
            self.x0[ix_['sigr4']] = float(2.8288294334)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr4']),float(2.8288294334)))
        if('x4' in ix_):
            self.x0[ix_['x4']] = float(60.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x4']),float(60.000000000)))
        if('y4' in ix_):
            self.x0[ix_['y4']] = float(4109.9328100)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y4']),float(4109.9328100)))
        if('w5' in ix_):
            self.x0[ix_['w5']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w5']),float(30.000000000)))
        if('sigt5' in ix_):
            self.x0[ix_['sigt5']] = float(66.059628152)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt5']),float(66.059628152)))
        if('sigr5' in ix_):
            self.x0[ix_['sigr5']] = float(3.4665315686)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr5']),float(3.4665315686)))
        if('x5' in ix_):
            self.x0[ix_['x5']] = float(75.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x5']),float(75.000000000)))
        if('y5' in ix_):
            self.x0[ix_['y5']] = float(5106.6886884)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y5']),float(5106.6886884)))
        if('w6' in ix_):
            self.x0[ix_['w6']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w6']),float(30.000000000)))
        if('sigt6' in ix_):
            self.x0[ix_['sigt6']] = float(65.301335165)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt6']),float(65.301335165)))
        if('sigr6' in ix_):
            self.x0[ix_['sigr6']] = float(4.0788400842)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr6']),float(4.0788400842)))
        if('x6' in ix_):
            self.x0[ix_['x6']] = float(90.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x6']),float(90.000000000)))
        if('y6' in ix_):
            self.x0[ix_['y6']] = float(6091.8959133)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y6']),float(6091.8959133)))
        if('w7' in ix_):
            self.x0[ix_['w7']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w7']),float(30.000000000)))
        if('sigt7' in ix_):
            self.x0[ix_['sigt7']] = float(64.565172618)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt7']),float(64.565172618)))
        if('sigr7' in ix_):
            self.x0[ix_['sigr7']] = float(4.6668413532)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr7']),float(4.6668413532)))
        if('x7' in ix_):
            self.x0[ix_['x7']] = float(105.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x7']),float(105.00000000)))
        if('y7' in ix_):
            self.x0[ix_['y7']] = float(7065.8947217)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y7']),float(7065.8947217)))
        if('w8' in ix_):
            self.x0[ix_['w8']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w8']),float(30.000000000)))
        if('sigt8' in ix_):
            self.x0[ix_['sigt8']] = float(63.850096500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt8']),float(63.850096500)))
        if('sigr8' in ix_):
            self.x0[ix_['sigr8']] = float(5.2315615316)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr8']),float(5.2315615316)))
        if('x8' in ix_):
            self.x0[ix_['x8']] = float(120.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x8']),float(120.00000000)))
        if('y8' in ix_):
            self.x0[ix_['y8']] = float(8029.0092401)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y8']),float(8029.0092401)))
        if('w9' in ix_):
            self.x0[ix_['w9']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w9']),float(30.000000000)))
        if('sigt9' in ix_):
            self.x0[ix_['sigt9']] = float(63.155118892)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt9']),float(63.155118892)))
        if('sigr9' in ix_):
            self.x0[ix_['sigr9']] = float(5.7739705049)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr9']),float(5.7739705049)))
        if('x9' in ix_):
            self.x0[ix_['x9']] = float(135.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x9']),float(135.00000000)))
        if('y9' in ix_):
            self.x0[ix_['y9']] = float(8981.5483555)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y9']),float(8981.5483555)))
        if('w10' in ix_):
            self.x0[ix_['w10']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w10']),float(30.000000000)))
        if('sigt10' in ix_):
            self.x0[ix_['sigt10']] = float(62.479304325)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt10']),float(62.479304325)))
        if('sigr10' in ix_):
            self.x0[ix_['sigr10']] = float(6.2949855370)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr10']),float(6.2949855370)))
        if('x10' in ix_):
            self.x0[ix_['x10']] = float(150.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x10']),float(150.00000000)))
        if('y10' in ix_):
            self.x0[ix_['y10']] = float(9923.8065296)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y10']),float(9923.8065296)))
        if('w11' in ix_):
            self.x0[ix_['w11']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w11']),float(30.000000000)))
        if('sigt11' in ix_):
            self.x0[ix_['sigt11']] = float(61.821766398)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt11']),float(61.821766398)))
        if('sigr11' in ix_):
            self.x0[ix_['sigr11']] = float(6.7954746441)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr11']),float(6.7954746441)))
        if('x11' in ix_):
            self.x0[ix_['x11']] = float(165.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x11']),float(165.00000000)))
        if('y11' in ix_):
            self.x0[ix_['y11']] = float(10856.064560)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y11']),float(10856.064560)))
        if('w12' in ix_):
            self.x0[ix_['w12']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w12']),float(30.000000000)))
        if('sigt12' in ix_):
            self.x0[ix_['sigt12']] = float(61.181664651)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt12']),float(61.181664651)))
        if('sigr12' in ix_):
            self.x0[ix_['sigr12']] = float(7.2762597204)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr12']),float(7.2762597204)))
        if('x12' in ix_):
            self.x0[ix_['x12']] = float(180.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x12']),float(180.00000000)))
        if('y12' in ix_):
            self.x0[ix_['y12']] = float(11778.590293)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y12']),float(11778.590293)))
        if('w13' in ix_):
            self.x0[ix_['w13']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w13']),float(30.000000000)))
        if('sigt13' in ix_):
            self.x0[ix_['sigt13']] = float(60.558201672)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt13']),float(60.558201672)))
        if('sigr13' in ix_):
            self.x0[ix_['sigr13']] = float(7.7381194336)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr13']),float(7.7381194336)))
        if('x13' in ix_):
            self.x0[ix_['x13']] = float(195.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x13']),float(195.00000000)))
        if('y13' in ix_):
            self.x0[ix_['y13']] = float(12691.639290)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y13']),float(12691.639290)))
        if('w14' in ix_):
            self.x0[ix_['w14']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w14']),float(30.000000000)))
        if('sigt14' in ix_):
            self.x0[ix_['sigt14']] = float(59.950620407)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt14']),float(59.950620407)))
        if('sigr14' in ix_):
            self.x0[ix_['sigr14']] = float(8.1817919107)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr14']),float(8.1817919107)))
        if('x14' in ix_):
            self.x0[ix_['x14']] = float(210.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x14']),float(210.00000000)))
        if('y14' in ix_):
            self.x0[ix_['y14']] = float(13595.455456)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y14']),float(13595.455456)))
        if('w15' in ix_):
            self.x0[ix_['w15']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w15']),float(30.000000000)))
        if('sigt15' in ix_):
            self.x0[ix_['sigt15']] = float(59.358201664)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt15']),float(59.358201664)))
        if('sigr15' in ix_):
            self.x0[ix_['sigr15']] = float(8.6079772309)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr15']),float(8.6079772309)))
        if('x15' in ix_):
            self.x0[ix_['x15']] = float(225.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x15']),float(225.00000000)))
        if('y15' in ix_):
            self.x0[ix_['y15']] = float(14490.271621)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y15']),float(14490.271621)))
        if('w16' in ix_):
            self.x0[ix_['w16']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w16']),float(30.000000000)))
        if('sigt16' in ix_):
            self.x0[ix_['sigt16']] = float(58.780261800)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt16']),float(58.780261800)))
        if('sigr16' in ix_):
            self.x0[ix_['sigr16']] = float(9.0173397415)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr16']),float(9.0173397415)))
        if('x16' in ix_):
            self.x0[ix_['x16']] = float(240.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x16']),float(240.00000000)))
        if('y16' in ix_):
            self.x0[ix_['y16']] = float(15376.310097)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y16']),float(15376.310097)))
        if('w17' in ix_):
            self.x0[ix_['w17']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w17']),float(30.000000000)))
        if('sigt17' in ix_):
            self.x0[ix_['sigt17']] = float(58.216150564)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt17']),float(58.216150564)))
        if('sigr17' in ix_):
            self.x0[ix_['sigr17']] = float(9.4105102106)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr17']),float(9.4105102106)))
        if('x17' in ix_):
            self.x0[ix_['x17']] = float(255.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x17']),float(255.00000000)))
        if('y17' in ix_):
            self.x0[ix_['y17']] = float(16253.783190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y17']),float(16253.783190)))
        if('w18' in ix_):
            self.x0[ix_['w18']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w18']),float(30.000000000)))
        if('sigt18' in ix_):
            self.x0[ix_['sigt18']] = float(57.665249095)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt18']),float(57.665249095)))
        if('sigr18' in ix_):
            self.x0[ix_['sigr18']] = float(9.7880878301)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr18']),float(9.7880878301)))
        if('x18' in ix_):
            self.x0[ix_['x18']] = float(270.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x18']),float(270.00000000)))
        if('y18' in ix_):
            self.x0[ix_['y18']] = float(17122.893688)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y18']),float(17122.893688)))
        if('w19' in ix_):
            self.x0[ix_['w19']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w19']),float(30.000000000)))
        if('sigt19' in ix_):
            self.x0[ix_['sigt19']] = float(57.126968056)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt19']),float(57.126968056)))
        if('sigr19' in ix_):
            self.x0[ix_['sigr19']] = float(10.150642080)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr19']),float(10.150642080)))
        if('x19' in ix_):
            self.x0[ix_['x19']] = float(285.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x19']),float(285.00000000)))
        if('y19' in ix_):
            self.x0[ix_['y19']] = float(17983.835316)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y19']),float(17983.835316)))
        if('w20' in ix_):
            self.x0[ix_['w20']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w20']),float(30.000000000)))
        if('sigt20' in ix_):
            self.x0[ix_['sigt20']] = float(56.600745894)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt20']),float(56.600745894)))
        if('sigr20' in ix_):
            self.x0[ix_['sigr20']] = float(10.498714467)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr20']),float(10.498714467)))
        if('x20' in ix_):
            self.x0[ix_['x20']] = float(300.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x20']),float(300.00000000)))
        if('y20' in ix_):
            self.x0[ix_['y20']] = float(18836.793171)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y20']),float(18836.793171)))
        if('w21' in ix_):
            self.x0[ix_['w21']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w21']),float(30.000000000)))
        if('sigt21' in ix_):
            self.x0[ix_['sigt21']] = float(56.086047224)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt21']),float(56.086047224)))
        if('sigr21' in ix_):
            self.x0[ix_['sigr21']] = float(10.832820143)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr21']),float(10.832820143)))
        if('x21' in ix_):
            self.x0[ix_['x21']] = float(315.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x21']),float(315.00000000)))
        if('y21' in ix_):
            self.x0[ix_['y21']] = float(19681.944119)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y21']),float(19681.944119)))
        if('w22' in ix_):
            self.x0[ix_['w22']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w22']),float(30.000000000)))
        if('sigt22' in ix_):
            self.x0[ix_['sigt22']] = float(55.582361311)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt22']),float(55.582361311)))
        if('sigr22' in ix_):
            self.x0[ix_['sigr22']] = float(11.153449416)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr22']),float(11.153449416)))
        if('x22' in ix_):
            self.x0[ix_['x22']] = float(330.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x22']),float(330.00000000)))
        if('y22' in ix_):
            self.x0[ix_['y22']] = float(20519.457183)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y22']),float(20519.457183)))
        if('w23' in ix_):
            self.x0[ix_['w23']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w23']),float(30.000000000)))
        if('sigt23' in ix_):
            self.x0[ix_['sigt23']] = float(55.089200663)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt23']),float(55.089200663)))
        if('sigr23' in ix_):
            self.x0[ix_['sigr23']] = float(11.461069163)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr23']),float(11.461069163)))
        if('x23' in ix_):
            self.x0[ix_['x23']] = float(345.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x23']),float(345.00000000)))
        if('y23' in ix_):
            self.x0[ix_['y23']] = float(21349.493898)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y23']),float(21349.493898)))
        if('w24' in ix_):
            self.x0[ix_['w24']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w24']),float(30.000000000)))
        if('sigt24' in ix_):
            self.x0[ix_['sigt24']] = float(54.606099709)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt24']),float(54.606099709)))
        if('sigr24' in ix_):
            self.x0[ix_['sigr24']] = float(11.756124148)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr24']),float(11.756124148)))
        if('x24' in ix_):
            self.x0[ix_['x24']] = float(360.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x24']),float(360.00000000)))
        if('y24' in ix_):
            self.x0[ix_['y24']] = float(22172.208651)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y24']),float(22172.208651)))
        if('w25' in ix_):
            self.x0[ix_['w25']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w25']),float(30.000000000)))
        if('sigt25' in ix_):
            self.x0[ix_['sigt25']] = float(54.132613569)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt25']),float(54.132613569)))
        if('sigr25' in ix_):
            self.x0[ix_['sigr25']] = float(12.039038253)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr25']),float(12.039038253)))
        if('x25' in ix_):
            self.x0[ix_['x25']] = float(375.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x25']),float(375.00000000)))
        if('y25' in ix_):
            self.x0[ix_['y25']] = float(22987.749000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y25']),float(22987.749000)))
        if('w26' in ix_):
            self.x0[ix_['w26']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w26']),float(30.000000000)))
        if('sigt26' in ix_):
            self.x0[ix_['sigt26']] = float(53.668316898)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt26']),float(53.668316898)))
        if('sigr26' in ix_):
            self.x0[ix_['sigr26']] = float(12.310215632)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr26']),float(12.310215632)))
        if('x26' in ix_):
            self.x0[ix_['x26']] = float(390.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x26']),float(390.00000000)))
        if('y26' in ix_):
            self.x0[ix_['y26']] = float(23796.255979)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y26']),float(23796.255979)))
        if('w27' in ix_):
            self.x0[ix_['w27']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w27']),float(30.000000000)))
        if('sigt27' in ix_):
            self.x0[ix_['sigt27']] = float(53.212802809)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt27']),float(53.212802809)))
        if('sigr27' in ix_):
            self.x0[ix_['sigr27']] = float(12.570041791)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr27']),float(12.570041791)))
        if('x27' in ix_):
            self.x0[ix_['x27']] = float(405.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x27']),float(405.00000000)))
        if('y27' in ix_):
            self.x0[ix_['y27']] = float(24597.864377)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y27']),float(24597.864377)))
        if('w28' in ix_):
            self.x0[ix_['w28']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w28']),float(30.000000000)))
        if('sigt28' in ix_):
            self.x0[ix_['sigt28']] = float(52.765681856)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt28']),float(52.765681856)))
        if('sigr28' in ix_):
            self.x0[ix_['sigr28']] = float(12.818884596)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr28']),float(12.818884596)))
        if('x28' in ix_):
            self.x0[ix_['x28']] = float(420.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x28']),float(420.00000000)))
        if('y28' in ix_):
            self.x0[ix_['y28']] = float(25392.703012)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y28']),float(25392.703012)))
        if('w29' in ix_):
            self.x0[ix_['w29']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w29']),float(30.000000000)))
        if('sigt29' in ix_):
            self.x0[ix_['sigt29']] = float(52.326581096)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt29']),float(52.326581096)))
        if('sigr29' in ix_):
            self.x0[ix_['sigr29']] = float(13.057095224)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr29']),float(13.057095224)))
        if('x29' in ix_):
            self.x0[ix_['x29']] = float(435.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x29']),float(435.00000000)))
        if('y29' in ix_):
            self.x0[ix_['y29']] = float(26180.894984)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y29']),float(26180.894984)))
        if('w30' in ix_):
            self.x0[ix_['w30']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w30']),float(30.000000000)))
        if('sigt30' in ix_):
            self.x0[ix_['sigt30']] = float(51.895143190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt30']),float(51.895143190)))
        if('sigr30' in ix_):
            self.x0[ix_['sigr30']] = float(13.285009049)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr30']),float(13.285009049)))
        if('x30' in ix_):
            self.x0[ix_['x30']] = float(450.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x30']),float(450.00000000)))
        if('y30' in ix_):
            self.x0[ix_['y30']] = float(26962.557916)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y30']),float(26962.557916)))
        if('w31' in ix_):
            self.x0[ix_['w31']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w31']),float(30.000000000)))
        if('sigt31' in ix_):
            self.x0[ix_['sigt31']] = float(51.471025575)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt31']),float(51.471025575)))
        if('sigr31' in ix_):
            self.x0[ix_['sigr31']] = float(13.502946477)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr31']),float(13.502946477)))
        if('x31' in ix_):
            self.x0[ix_['x31']] = float(465.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x31']),float(465.00000000)))
        if('y31' in ix_):
            self.x0[ix_['y31']] = float(27737.804182)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y31']),float(27737.804182)))
        if('w32' in ix_):
            self.x0[ix_['w32']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w32']),float(30.000000000)))
        if('sigt32' in ix_):
            self.x0[ix_['sigt32']] = float(51.053899680)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt32']),float(51.053899680)))
        if('sigr32' in ix_):
            self.x0[ix_['sigr32']] = float(13.711213727)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr32']),float(13.711213727)))
        if('x32' in ix_):
            self.x0[ix_['x32']] = float(480.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x32']),float(480.00000000)))
        if('y32' in ix_):
            self.x0[ix_['y32']] = float(28506.741121)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y32']),float(28506.741121)))
        if('w33' in ix_):
            self.x0[ix_['w33']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w33']),float(30.000000000)))
        if('sigt33' in ix_):
            self.x0[ix_['sigt33']] = float(50.643450190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt33']),float(50.643450190)))
        if('sigr33' in ix_):
            self.x0[ix_['sigr33']] = float(13.910103570)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr33']),float(13.910103570)))
        if('x33' in ix_):
            self.x0[ix_['x33']] = float(495.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x33']),float(495.00000000)))
        if('y33' in ix_):
            self.x0[ix_['y33']] = float(29269.471245)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y33']),float(29269.471245)))
        if('w34' in ix_):
            self.x0[ix_['w34']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w34']),float(30.000000000)))
        if('sigt34' in ix_):
            self.x0[ix_['sigt34']] = float(50.239374354)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt34']),float(50.239374354)))
        if('sigr34' in ix_):
            self.x0[ix_['sigr34']] = float(14.099896013)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr34']),float(14.099896013)))
        if('x34' in ix_):
            self.x0[ix_['x34']] = float(510.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x34']),float(510.00000000)))
        if('y34' in ix_):
            self.x0[ix_['y34']] = float(30026.092429)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y34']),float(30026.092429)))
        if('w35' in ix_):
            self.x0[ix_['w35']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w35']),float(30.000000000)))
        if('sigt35' in ix_):
            self.x0[ix_['sigt35']] = float(49.841381338)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt35']),float(49.841381338)))
        if('sigr35' in ix_):
            self.x0[ix_['sigr35']] = float(14.280858959)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr35']),float(14.280858959)))
        if('x35' in ix_):
            self.x0[ix_['x35']] = float(525.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x35']),float(525.00000000)))
        if('y35' in ix_):
            self.x0[ix_['y35']] = float(30776.698097)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y35']),float(30776.698097)))
        if('w36' in ix_):
            self.x0[ix_['w36']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w36']),float(30.000000000)))
        if('sigt36' in ix_):
            self.x0[ix_['sigt36']] = float(49.449191607)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt36']),float(49.449191607)))
        if('sigr36' in ix_):
            self.x0[ix_['sigr36']] = float(14.453248807)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr36']),float(14.453248807)))
        if('x36' in ix_):
            self.x0[ix_['x36']] = float(540.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x36']),float(540.00000000)))
        if('y36' in ix_):
            self.x0[ix_['y36']] = float(31521.377394)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y36']),float(31521.377394)))
        if('w37' in ix_):
            self.x0[ix_['w37']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w37']),float(30.000000000)))
        if('sigt37' in ix_):
            self.x0[ix_['sigt37']] = float(49.062536356)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt37']),float(49.062536356)))
        if('sigr37' in ix_):
            self.x0[ix_['sigr37']] = float(14.617311038)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr37']),float(14.617311038)))
        if('x37' in ix_):
            self.x0[ix_['x37']] = float(555.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x37']),float(555.00000000)))
        if('y37' in ix_):
            self.x0[ix_['y37']] = float(32260.215354)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y37']),float(32260.215354)))
        if('w38' in ix_):
            self.x0[ix_['w38']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w38']),float(30.000000000)))
        if('sigt38' in ix_):
            self.x0[ix_['sigt38']] = float(48.681156965)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt38']),float(48.681156965)))
        if('sigr38' in ix_):
            self.x0[ix_['sigr38']] = float(14.773280750)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr38']),float(14.773280750)))
        if('x38' in ix_):
            self.x0[ix_['x38']] = float(570.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x38']),float(570.00000000)))
        if('y38' in ix_):
            self.x0[ix_['y38']] = float(32993.293054)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y38']),float(32993.293054)))
        if('w39' in ix_):
            self.x0[ix_['w39']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w39']),float(30.000000000)))
        if('sigt39' in ix_):
            self.x0[ix_['sigt39']] = float(48.304804485)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt39']),float(48.304804485)))
        if('sigr39' in ix_):
            self.x0[ix_['sigr39']] = float(14.921383174)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr39']),float(14.921383174)))
        if('x39' in ix_):
            self.x0[ix_['x39']] = float(585.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x39']),float(585.00000000)))
        if('y39' in ix_):
            self.x0[ix_['y39']] = float(33720.687765)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y39']),float(33720.687765)))
        if('w40' in ix_):
            self.x0[ix_['w40']] = float(30.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w40']),float(30.000000000)))
        if('sigt40' in ix_):
            self.x0[ix_['sigt40']] = float(47.933239160)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt40']),float(47.933239160)))
        if('sigr40' in ix_):
            self.x0[ix_['sigr40']] = float(15.061834152)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr40']),float(15.061834152)))
        if('x40' in ix_):
            self.x0[ix_['x40']] = float(600.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x40']),float(600.00000000)))
        if('y40' in ix_):
            self.x0[ix_['y40']] = float(34442.473092)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y40']),float(34442.473092)))
        if('w41' in ix_):
            self.x0[ix_['w41']] = float(29.017500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w41']),float(29.017500000)))
        if('sigt41' in ix_):
            self.x0[ix_['sigt41']] = float(47.721358592)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt41']),float(47.721358592)))
        if('sigr41' in ix_):
            self.x0[ix_['sigr41']] = float(15.705670256)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr41']),float(15.705670256)))
        if('x41' in ix_):
            self.x0[ix_['x41']] = float(614.75437500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x41']),float(614.75437500)))
        if('y41' in ix_):
            self.x0[ix_['y41']] = float(35148.161016)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y41']),float(35148.161016)))
        if('w42' in ix_):
            self.x0[ix_['w42']] = float(28.035000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w42']),float(28.035000000)))
        if('sigt42' in ix_):
            self.x0[ix_['sigt42']] = float(47.528871380)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt42']),float(47.528871380)))
        if('sigr42' in ix_):
            self.x0[ix_['sigr42']] = float(16.379639595)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr42']),float(16.379639595)))
        if('x42' in ix_):
            self.x0[ix_['x42']] = float(629.01750000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x42']),float(629.01750000)))
        if('y42' in ix_):
            self.x0[ix_['y42']] = float(35827.467624)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y42']),float(35827.467624)))
        if('w43' in ix_):
            self.x0[ix_['w43']] = float(27.052500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w43']),float(27.052500000)))
        if('sigt43' in ix_):
            self.x0[ix_['sigt43']] = float(47.356912245)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt43']),float(47.356912245)))
        if('sigr43' in ix_):
            self.x0[ix_['sigr43']] = float(17.087816055)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr43']),float(17.087816055)))
        if('x43' in ix_):
            self.x0[ix_['x43']] = float(642.78937500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x43']),float(642.78937500)))
        if('y43' in ix_):
            self.x0[ix_['y43']] = float(36480.866319)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y43']),float(36480.866319)))
        if('w44' in ix_):
            self.x0[ix_['w44']] = float(26.070000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w44']),float(26.070000000)))
        if('sigt44' in ix_):
            self.x0[ix_['sigt44']] = float(47.206824226)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt44']),float(47.206824226)))
        if('sigr44' in ix_):
            self.x0[ix_['sigr44']] = float(17.834855648)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr44']),float(17.834855648)))
        if('x44' in ix_):
            self.x0[ix_['x44']] = float(656.07000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x44']),float(656.07000000)))
        if('y44' in ix_):
            self.x0[ix_['y44']] = float(37108.817513)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y44']),float(37108.817513)))
        if('w45' in ix_):
            self.x0[ix_['w45']] = float(25.087500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w45']),float(25.087500000)))
        if('sigt45' in ix_):
            self.x0[ix_['sigt45']] = float(47.080196532)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt45']),float(47.080196532)))
        if('sigr45' in ix_):
            self.x0[ix_['sigr45']] = float(18.626112102)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr45']),float(18.626112102)))
        if('x45' in ix_):
            self.x0[ix_['x45']] = float(668.85937500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x45']),float(668.85937500)))
        if('y45' in ix_):
            self.x0[ix_['y45']] = float(37711.769097)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y45']),float(37711.769097)))
        if('w46' in ix_):
            self.x0[ix_['w46']] = float(24.105000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w46']),float(24.105000000)))
        if('sigt46' in ix_):
            self.x0[ix_['sigt46']] = float(46.978911596)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt46']),float(46.978911596)))
        if('sigr46' in ix_):
            self.x0[ix_['sigr46']] = float(19.467780620)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr46']),float(19.467780620)))
        if('x46' in ix_):
            self.x0[ix_['x46']] = float(681.15750000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x46']),float(681.15750000)))
        if('y46' in ix_):
            self.x0[ix_['y46']] = float(38290.156871)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y46']),float(38290.156871)))
        if('w47' in ix_):
            self.x0[ix_['w47']] = float(23.122500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w47']),float(23.122500000)))
        if('sigt47' in ix_):
            self.x0[ix_['sigt47']] = float(46.905204055)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt47']),float(46.905204055)))
        if('sigr47' in ix_):
            self.x0[ix_['sigr47']] = float(20.367078186)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr47']),float(20.367078186)))
        if('x47' in ix_):
            self.x0[ix_['x47']] = float(692.96437500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x47']),float(692.96437500)))
        if('y47' in ix_):
            self.x0[ix_['y47']] = float(38844.404932)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y47']),float(38844.404932)))
        if('w48' in ix_):
            self.x0[ix_['w48']] = float(22.140000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w48']),float(22.140000000)))
        if('sigt48' in ix_):
            self.x0[ix_['sigt48']] = float(46.861735282)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt48']),float(46.861735282)))
        if('sigr48' in ix_):
            self.x0[ix_['sigr48']] = float(21.332471779)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr48']),float(21.332471779)))
        if('x48' in ix_):
            self.x0[ix_['x48']] = float(704.28000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x48']),float(704.28000000)))
        if('y48' in ix_):
            self.x0[ix_['y48']] = float(39374.926032)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y48']),float(39374.926032)))
        if('w49' in ix_):
            self.x0[ix_['w49']] = float(21.372500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w49']),float(21.372500000)))
        if('sigt49' in ix_):
            self.x0[ix_['sigt49']] = float(46.783637776)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt49']),float(46.783637776)))
        if('sigr49' in ix_):
            self.x0[ix_['sigr49']] = float(22.149717856)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr49']),float(22.149717856)))
        if('x49' in ix_):
            self.x0[ix_['x49']] = float(715.15812500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x49']),float(715.15812500)))
        if('y49' in ix_):
            self.x0[ix_['y49']] = float(39884.276561)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y49']),float(39884.276561)))
        if('w50' in ix_):
            self.x0[ix_['w50']] = float(20.605000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w50']),float(20.605000000)))
        if('sigt50' in ix_):
            self.x0[ix_['sigt50']] = float(46.729531640)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt50']),float(46.729531640)))
        if('sigr50' in ix_):
            self.x0[ix_['sigr50']] = float(23.016346353)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr50']),float(23.016346353)))
        if('x50' in ix_):
            self.x0[ix_['x50']] = float(725.65250000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x50']),float(725.65250000)))
        if('y50' in ix_):
            self.x0[ix_['y50']] = float(40374.962886)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y50']),float(40374.962886)))
        if('w51' in ix_):
            self.x0[ix_['w51']] = float(19.837500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w51']),float(19.837500000)))
        if('sigt51' in ix_):
            self.x0[ix_['sigt51']] = float(46.701411241)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt51']),float(46.701411241)))
        if('sigr51' in ix_):
            self.x0[ix_['sigr51']] = float(23.938737423)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr51']),float(23.938737423)))
        if('x51' in ix_):
            self.x0[ix_['x51']] = float(735.76312500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x51']),float(735.76312500)))
        if('y51' in ix_):
            self.x0[ix_['y51']] = float(40847.288197)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y51']),float(40847.288197)))
        if('w52' in ix_):
            self.x0[ix_['w52']] = float(19.070000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w52']),float(19.070000000)))
        if('sigt52' in ix_):
            self.x0[ix_['sigt52']] = float(46.701615132)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt52']),float(46.701615132)))
        if('sigr52' in ix_):
            self.x0[ix_['sigr52']] = float(24.924274110)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr52']),float(24.924274110)))
        if('x52' in ix_):
            self.x0[ix_['x52']] = float(745.49000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x52']),float(745.49000000)))
        if('y52' in ix_):
            self.x0[ix_['y52']] = float(41301.547959)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y52']),float(41301.547959)))
        if('w53' in ix_):
            self.x0[ix_['w53']] = float(18.302500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w53']),float(18.302500000)))
        if('sigt53' in ix_):
            self.x0[ix_['sigt53']] = float(46.732895222)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt53']),float(46.732895222)))
        if('sigr53' in ix_):
            self.x0[ix_['sigr53']] = float(25.981553727)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr53']),float(25.981553727)))
        if('x53' in ix_):
            self.x0[ix_['x53']] = float(754.83312500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x53']),float(754.83312500)))
        if('y53' in ix_):
            self.x0[ix_['y53']] = float(41738.030113)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y53']),float(41738.030113)))
        if('w54' in ix_):
            self.x0[ix_['w54']] = float(17.535000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w54']),float(17.535000000)))
        if('sigt54' in ix_):
            self.x0[ix_['sigt54']] = float(46.798503908)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt54']),float(46.798503908)))
        if('sigr54' in ix_):
            self.x0[ix_['sigr54']] = float(27.120654675)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr54']),float(27.120654675)))
        if('x54' in ix_):
            self.x0[ix_['x54']] = float(763.79250000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x54']),float(763.79250000)))
        if('y54' in ix_):
            self.x0[ix_['y54']] = float(42157.015258)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y54']),float(42157.015258)))
        if('w55' in ix_):
            self.x0[ix_['w55']] = float(16.767500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w55']),float(16.767500000)))
        if('sigt55' in ix_):
            self.x0[ix_['sigt55']] = float(46.902304859)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt55']),float(46.902304859)))
        if('sigr55' in ix_):
            self.x0[ix_['sigr55']] = float(28.353476460)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr55']),float(28.353476460)))
        if('x55' in ix_):
            self.x0[ix_['x55']] = float(772.36812500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x55']),float(772.36812500)))
        if('y55' in ix_):
            self.x0[ix_['y55']] = float(42558.776798)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y55']),float(42558.776798)))
        if('w56' in ix_):
            self.x0[ix_['w56']] = float(16.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w56']),float(16.000000000)))
        if('sigt56' in ix_):
            self.x0[ix_['sigt56']] = float(47.048915288)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt56']),float(47.048915288)))
        if('sigr56' in ix_):
            self.x0[ix_['sigr56']] = float(29.694177498)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr56']),float(29.694177498)))
        if('x56' in ix_):
            self.x0[ix_['x56']] = float(780.56000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x56']),float(780.56000000)))
        if('y56' in ix_):
            self.x0[ix_['y56']] = float(42943.581059)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y56']),float(42943.581059)))
        if('w57' in ix_):
            self.x0[ix_['w57']] = float(15.440000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w57']),float(15.440000000)))
        if('sigt57' in ix_):
            self.x0[ix_['sigt57']] = float(47.117144627)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt57']),float(47.117144627)))
        if('sigr57' in ix_):
            self.x0[ix_['sigr57']] = float(30.741797333)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr57']),float(30.741797333)))
        if('x57' in ix_):
            self.x0[ix_['x57']] = float(788.42000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x57']),float(788.42000000)))
        if('y57' in ix_):
            self.x0[ix_['y57']] = float(43313.648898)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y57']),float(43313.648898)))
        if('w58' in ix_):
            self.x0[ix_['w58']] = float(14.880000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w58']),float(14.880000000)))
        if('sigt58' in ix_):
            self.x0[ix_['sigt58']] = float(47.215105046)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt58']),float(47.215105046)))
        if('sigr58' in ix_):
            self.x0[ix_['sigr58']] = float(31.860061245)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr58']),float(31.860061245)))
        if('x58' in ix_):
            self.x0[ix_['x58']] = float(796.00000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x58']),float(796.00000000)))
        if('y58' in ix_):
            self.x0[ix_['y58']] = float(43671.161267)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y58']),float(43671.161267)))
        if('w59' in ix_):
            self.x0[ix_['w59']] = float(14.320000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w59']),float(14.320000000)))
        if('sigt59' in ix_):
            self.x0[ix_['sigt59']] = float(47.345666162)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt59']),float(47.345666162)))
        if('sigr59' in ix_):
            self.x0[ix_['sigr59']] = float(33.057769694)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr59']),float(33.057769694)))
        if('x59' in ix_):
            self.x0[ix_['x59']] = float(803.30000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x59']),float(803.30000000)))
        if('y59' in ix_):
            self.x0[ix_['y59']] = float(44016.298943)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y59']),float(44016.298943)))
        if('w60' in ix_):
            self.x0[ix_['w60']] = float(13.760000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w60']),float(13.760000000)))
        if('sigt60' in ix_):
            self.x0[ix_['sigt60']] = float(47.512175218)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt60']),float(47.512175218)))
        if('sigr60' in ix_):
            self.x0[ix_['sigr60']] = float(34.345138122)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr60']),float(34.345138122)))
        if('x60' in ix_):
            self.x0[ix_['x60']] = float(810.32000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x60']),float(810.32000000)))
        if('y60' in ix_):
            self.x0[ix_['y60']] = float(44349.238310)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y60']),float(44349.238310)))
        if('w61' in ix_):
            self.x0[ix_['w61']] = float(13.200000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w61']),float(13.200000000)))
        if('sigt61' in ix_):
            self.x0[ix_['sigt61']] = float(47.718555318)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt61']),float(47.718555318)))
        if('sigr61' in ix_):
            self.x0[ix_['sigr61']] = float(35.734097806)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr61']),float(35.734097806)))
        if('x61' in ix_):
            self.x0[ix_['x61']] = float(817.06000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x61']),float(817.06000000)))
        if('y61' in ix_):
            self.x0[ix_['y61']] = float(44670.151426)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y61']),float(44670.151426)))
        if('w62' in ix_):
            self.x0[ix_['w62']] = float(12.640000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w62']),float(12.640000000)))
        if('sigt62' in ix_):
            self.x0[ix_['sigt62']] = float(47.969429433)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt62']),float(47.969429433)))
        if('sigr62' in ix_):
            self.x0[ix_['sigr62']] = float(37.238676641)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr62']),float(37.238676641)))
        if('x62' in ix_):
            self.x0[ix_['x62']] = float(823.52000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x62']),float(823.52000000)))
        if('y62' in ix_):
            self.x0[ix_['y62']] = float(44979.206055)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y62']),float(44979.206055)))
        if('w63' in ix_):
            self.x0[ix_['w63']] = float(12.080000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w63']),float(12.080000000)))
        if('sigt63' in ix_):
            self.x0[ix_['sigt63']] = float(48.270278469)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt63']),float(48.270278469)))
        if('sigr63' in ix_):
            self.x0[ix_['sigr63']] = float(38.875485767)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr63']),float(38.875485767)))
        if('x63' in ix_):
            self.x0[ix_['x63']] = float(829.70000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x63']),float(829.70000000)))
        if('y63' in ix_):
            self.x0[ix_['y63']] = float(45276.565693)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y63']),float(45276.565693)))
        if('w64' in ix_):
            self.x0[ix_['w64']] = float(11.520000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w64']),float(11.520000000)))
        if('sigt64' in ix_):
            self.x0[ix_['sigt64']] = float(48.627644833)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt64']),float(48.627644833)))
        if('sigr64' in ix_):
            self.x0[ix_['sigr64']] = float(40.664348076)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr64']),float(40.664348076)))
        if('x64' in ix_):
            self.x0[ix_['x64']] = float(835.60000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x64']),float(835.60000000)))
        if('y64' in ix_):
            self.x0[ix_['y64']] = float(45562.389551)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y64']),float(45562.389551)))
        if('w65' in ix_):
            self.x0[ix_['w65']] = float(11.175000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w65']),float(11.175000000)))
        if('sigt65' in ix_):
            self.x0[ix_['sigt65']] = float(48.801073466)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt65']),float(48.801073466)))
        if('sigr65' in ix_):
            self.x0[ix_['sigr65']] = float(41.809788300)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr65']),float(41.809788300)))
        if('x65' in ix_):
            self.x0[ix_['x65']] = float(841.27375000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x65']),float(841.27375000)))
        if('y65' in ix_):
            self.x0[ix_['y65']] = float(45838.775167)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y65']),float(45838.775167)))
        if('w66' in ix_):
            self.x0[ix_['w66']] = float(10.830000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w66']),float(10.830000000)))
        if('sigt66' in ix_):
            self.x0[ix_['sigt66']] = float(49.001998851)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt66']),float(49.001998851)))
        if('sigr66' in ix_):
            self.x0[ix_['sigr66']] = float(43.023368828)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr66']),float(43.023368828)))
        if('x66' in ix_):
            self.x0[ix_['x66']] = float(846.77500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x66']),float(846.77500000)))
        if('y66' in ix_):
            self.x0[ix_['y66']] = float(46107.786078)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y66']),float(46107.786078)))
        if('w67' in ix_):
            self.x0[ix_['w67']] = float(10.485000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w67']),float(10.485000000)))
        if('sigt67' in ix_):
            self.x0[ix_['sigt67']] = float(49.232764164)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt67']),float(49.232764164)))
        if('sigr67' in ix_):
            self.x0[ix_['sigr67']] = float(44.312155780)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr67']),float(44.312155780)))
        if('x67' in ix_):
            self.x0[ix_['x67']] = float(852.10375000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x67']),float(852.10375000)))
        if('y67' in ix_):
            self.x0[ix_['y67']] = float(46369.510373)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y67']),float(46369.510373)))
        if('w68' in ix_):
            self.x0[ix_['w68']] = float(10.140000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w68']),float(10.140000000)))
        if('sigt68' in ix_):
            self.x0[ix_['sigt68']] = float(49.496036475)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt68']),float(49.496036475)))
        if('sigr68' in ix_):
            self.x0[ix_['sigr68']] = float(45.684166505)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr68']),float(45.684166505)))
        if('x68' in ix_):
            self.x0[ix_['x68']] = float(857.26000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x68']),float(857.26000000)))
        if('y68' in ix_):
            self.x0[ix_['y68']] = float(46624.034209)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y68']),float(46624.034209)))
        if('w69' in ix_):
            self.x0[ix_['w69']] = float(9.7950000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w69']),float(9.7950000000)))
        if('sigt69' in ix_):
            self.x0[ix_['sigt69']] = float(49.794861993)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt69']),float(49.794861993)))
        if('sigr69' in ix_):
            self.x0[ix_['sigr69']] = float(47.148537492)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr69']),float(47.148537492)))
        if('x69' in ix_):
            self.x0[ix_['x69']] = float(862.24375000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x69']),float(862.24375000)))
        if('y69' in ix_):
            self.x0[ix_['y69']] = float(46871.441830)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y69']),float(46871.441830)))
        if('w70' in ix_):
            self.x0[ix_['w70']] = float(9.4500000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w70']),float(9.4500000000)))
        if('sigt70' in ix_):
            self.x0[ix_['sigt70']] = float(50.132733249)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt70']),float(50.132733249)))
        if('sigr70' in ix_):
            self.x0[ix_['sigr70']] = float(48.715729024)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr70']),float(48.715729024)))
        if('x70' in ix_):
            self.x0[ix_['x70']] = float(867.05500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x70']),float(867.05500000)))
        if('y70' in ix_):
            self.x0[ix_['y70']] = float(47111.815580)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y70']),float(47111.815580)))
        if('w71' in ix_):
            self.x0[ix_['w71']] = float(9.1050000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w71']),float(9.1050000000)))
        if('sigt71' in ix_):
            self.x0[ix_['sigt71']] = float(50.513671349)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt71']),float(50.513671349)))
        if('sigr71' in ix_):
            self.x0[ix_['sigr71']] = float(50.397776340)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr71']),float(50.397776340)))
        if('x71' in ix_):
            self.x0[ix_['x71']] = float(871.69375000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x71']),float(871.69375000)))
        if('y71' in ix_):
            self.x0[ix_['y71']] = float(47345.235907)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y71']),float(47345.235907)))
        if('w72' in ix_):
            self.x0[ix_['w72']] = float(8.7600000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w72']),float(8.7600000000)))
        if('sigt72' in ix_):
            self.x0[ix_['sigt72']] = float(50.942327392)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt72']),float(50.942327392)))
        if('sigr72' in ix_):
            self.x0[ix_['sigr72']] = float(52.208600106)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr72']),float(52.208600106)))
        if('x72' in ix_):
            self.x0[ix_['x72']] = float(876.16000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x72']),float(876.16000000)))
        if('y72' in ix_):
            self.x0[ix_['y72']] = float(47571.781348)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y72']),float(47571.781348)))
        if('w73' in ix_):
            self.x0[ix_['w73']] = float(8.6275000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w73']),float(8.6275000000)))
        if('sigt73' in ix_):
            self.x0[ix_['sigt73']] = float(51.020192955)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt73']),float(51.020192955)))
        if('sigr73' in ix_):
            self.x0[ix_['sigr73']] = float(52.831000764)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr73']),float(52.831000764)))
        if('x73' in ix_):
            self.x0[ix_['x73']] = float(880.50687500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x73']),float(880.50687500)))
        if('y73' in ix_):
            self.x0[ix_['y73']] = float(47793.389224)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y73']),float(47793.389224)))
        if('w74' in ix_):
            self.x0[ix_['w74']] = float(8.4950000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w74']),float(8.4950000000)))
        if('sigt74' in ix_):
            self.x0[ix_['sigt74']] = float(51.105469720)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt74']),float(51.105469720)))
        if('sigr74' in ix_):
            self.x0[ix_['sigr74']] = float(53.471000483)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr74']),float(53.471000483)))
        if('x74' in ix_):
            self.x0[ix_['x74']] = float(884.78750000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x74']),float(884.78750000)))
        if('y74' in ix_):
            self.x0[ix_['y74']] = float(48011.968644)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y74']),float(48011.968644)))
        if('w75' in ix_):
            self.x0[ix_['w75']] = float(8.3625000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w75']),float(8.3625000000)))
        if('sigt75' in ix_):
            self.x0[ix_['sigt75']] = float(51.198444013)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt75']),float(51.198444013)))
        if('sigr75' in ix_):
            self.x0[ix_['sigr75']] = float(54.129557211)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr75']),float(54.129557211)))
        if('x75' in ix_):
            self.x0[ix_['x75']] = float(889.00187500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x75']),float(889.00187500)))
        if('y75' in ix_):
            self.x0[ix_['y75']] = float(48227.540632)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y75']),float(48227.540632)))
        if('w76' in ix_):
            self.x0[ix_['w76']] = float(8.2300000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w76']),float(8.2300000000)))
        if('sigt76' in ix_):
            self.x0[ix_['sigt76']] = float(51.299424731)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt76']),float(51.299424731)))
        if('sigr76' in ix_):
            self.x0[ix_['sigr76']] = float(54.807687715)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr76']),float(54.807687715)))
        if('x76' in ix_):
            self.x0[ix_['x76']] = float(893.15000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x76']),float(893.15000000)))
        if('y76' in ix_):
            self.x0[ix_['y76']] = float(48440.125946)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y76']),float(48440.125946)))
        if('w77' in ix_):
            self.x0[ix_['w77']] = float(8.0975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w77']),float(8.0975000000)))
        if('sigt77' in ix_):
            self.x0[ix_['sigt77']] = float(51.408745006)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt77']),float(51.408745006)))
        if('sigr77' in ix_):
            self.x0[ix_['sigr77']] = float(55.506472516)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr77']),float(55.506472516)))
        if('x77' in ix_):
            self.x0[ix_['x77']] = float(897.23187500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x77']),float(897.23187500)))
        if('y77' in ix_):
            self.x0[ix_['y77']] = float(48649.745090)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y77']),float(48649.745090)))
        if('w78' in ix_):
            self.x0[ix_['w78']] = float(7.9650000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w78']),float(7.9650000000)))
        if('sigt78' in ix_):
            self.x0[ix_['sigt78']] = float(51.526764037)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt78']),float(51.526764037)))
        if('sigr78' in ix_):
            self.x0[ix_['sigr78']] = float(56.227061303)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr78']),float(56.227061303)))
        if('x78' in ix_):
            self.x0[ix_['x78']] = float(901.24750000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x78']),float(901.24750000)))
        if('y78' in ix_):
            self.x0[ix_['y78']] = float(48856.418337)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y78']),float(48856.418337)))
        if('w79' in ix_):
            self.x0[ix_['w79']] = float(7.8325000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w79']),float(7.8325000000)))
        if('sigt79' in ix_):
            self.x0[ix_['sigt79']] = float(51.653869098)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt79']),float(51.653869098)))
        if('sigr79' in ix_):
            self.x0[ix_['sigr79']] = float(56.970678892)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr79']),float(56.970678892)))
        if('x79' in ix_):
            self.x0[ix_['x79']] = float(905.19687500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x79']),float(905.19687500)))
        if('y79' in ix_):
            self.x0[ix_['y79']] = float(49060.165739)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y79']),float(49060.165739)))
        if('w80' in ix_):
            self.x0[ix_['w80']] = float(7.7000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w80']),float(7.7000000000)))
        if('sigt80' in ix_):
            self.x0[ix_['sigt80']] = float(51.790477770)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt80']),float(51.790477770)))
        if('sigr80' in ix_):
            self.x0[ix_['sigr80']] = float(57.738631802)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr80']),float(57.738631802)))
        if('x80' in ix_):
            self.x0[ix_['x80']] = float(909.08000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x80']),float(909.08000000)))
        if('y80' in ix_):
            self.x0[ix_['y80']] = float(49261.007141)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y80']),float(49261.007141)))
        if('w81' in ix_):
            self.x0[ix_['w81']] = float(7.6725000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w81']),float(7.6725000000)))
        if('sigt81' in ix_):
            self.x0[ix_['sigt81']] = float(51.694572152)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt81']),float(51.694572152)))
        if('sigr81' in ix_):
            self.x0[ix_['sigr81']] = float(57.731509685)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr81']),float(57.731509685)))
        if('x81' in ix_):
            self.x0[ix_['x81']] = float(912.92312500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x81']),float(912.92312500)))
        if('y81' in ix_):
            self.x0[ix_['y81']] = float(49459.860462)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y81']),float(49459.860462)))
        if('w82' in ix_):
            self.x0[ix_['w82']] = float(7.6450000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w82']),float(7.6450000000)))
        if('sigt82' in ix_):
            self.x0[ix_['sigt82']] = float(51.596245048)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt82']),float(51.596245048)))
        if('sigr82' in ix_):
            self.x0[ix_['sigr82']] = float(57.723691663)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr82']),float(57.723691663)))
        if('x82' in ix_):
            self.x0[ix_['x82']] = float(916.75250000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x82']),float(916.75250000)))
        if('y82' in ix_):
            self.x0[ix_['y82']] = float(49657.630436)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y82']),float(49657.630436)))
        if('w83' in ix_):
            self.x0[ix_['w83']] = float(7.6175000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w83']),float(7.6175000000)))
        if('sigt83' in ix_):
            self.x0[ix_['sigt83']] = float(51.495471481)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt83']),float(51.495471481)))
        if('sigr83' in ix_):
            self.x0[ix_['sigr83']] = float(57.715173744)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr83']),float(57.715173744)))
        if('x83' in ix_):
            self.x0[ix_['x83']] = float(920.56812500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x83']),float(920.56812500)))
        if('y83' in ix_):
            self.x0[ix_['y83']] = float(49854.310448)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y83']),float(49854.310448)))
        if('w84' in ix_):
            self.x0[ix_['w84']] = float(7.5900000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w84']),float(7.5900000000)))
        if('sigt84' in ix_):
            self.x0[ix_['sigt84']] = float(51.392226289)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt84']),float(51.392226289)))
        if('sigr84' in ix_):
            self.x0[ix_['sigr84']] = float(57.705951941)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr84']),float(57.705951941)))
        if('x84' in ix_):
            self.x0[ix_['x84']] = float(924.37000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x84']),float(924.37000000)))
        if('y84' in ix_):
            self.x0[ix_['y84']] = float(50049.893886)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y84']),float(50049.893886)))
        if('w85' in ix_):
            self.x0[ix_['w85']] = float(7.5625000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w85']),float(7.5625000000)))
        if('sigt85' in ix_):
            self.x0[ix_['sigt85']] = float(51.286484124)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt85']),float(51.286484124)))
        if('sigr85' in ix_):
            self.x0[ix_['sigr85']] = float(57.696022277)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr85']),float(57.696022277)))
        if('x85' in ix_):
            self.x0[ix_['x85']] = float(928.15812500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x85']),float(928.15812500)))
        if('y85' in ix_):
            self.x0[ix_['y85']] = float(50244.374144)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y85']),float(50244.374144)))
        if('w86' in ix_):
            self.x0[ix_['w86']] = float(7.5350000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w86']),float(7.5350000000)))
        if('sigt86' in ix_):
            self.x0[ix_['sigt86']] = float(51.178219453)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt86']),float(51.178219453)))
        if('sigr86' in ix_):
            self.x0[ix_['sigr86']] = float(57.685380778)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr86']),float(57.685380778)))
        if('x86' in ix_):
            self.x0[ix_['x86']] = float(931.93250000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x86']),float(931.93250000)))
        if('y86' in ix_):
            self.x0[ix_['y86']] = float(50437.744624)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y86']),float(50437.744624)))
        if('w87' in ix_):
            self.x0[ix_['w87']] = float(7.5075000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w87']),float(7.5075000000)))
        if('sigt87' in ix_):
            self.x0[ix_['sigt87']] = float(51.067406560)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt87']),float(51.067406560)))
        if('sigr87' in ix_):
            self.x0[ix_['sigr87']] = float(57.674023474)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr87']),float(57.674023474)))
        if('x87' in ix_):
            self.x0[ix_['x87']] = float(935.69312500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x87']),float(935.69312500)))
        if('y87' in ix_):
            self.x0[ix_['y87']] = float(50629.998734)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y87']),float(50629.998734)))
        if('w88' in ix_):
            self.x0[ix_['w88']] = float(7.4800000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w88']),float(7.4800000000)))
        if('sigt88' in ix_):
            self.x0[ix_['sigt88']] = float(50.954019546)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt88']),float(50.954019546)))
        if('sigr88' in ix_):
            self.x0[ix_['sigr88']] = float(57.661946402)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr88']),float(57.661946402)))
        if('x88' in ix_):
            self.x0[ix_['x88']] = float(939.44000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x88']),float(939.44000000)))
        if('y88' in ix_):
            self.x0[ix_['y88']] = float(50821.129889)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y88']),float(50821.129889)))
        if('w89' in ix_):
            self.x0[ix_['w89']] = float(7.4450000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w89']),float(7.4450000000)))
        if('sigt89' in ix_):
            self.x0[ix_['sigt89']] = float(50.855607376)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt89']),float(50.855607376)))
        if('sigr89' in ix_):
            self.x0[ix_['sigr89']] = float(57.707215988)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr89']),float(57.707215988)))
        if('x89' in ix_):
            self.x0[ix_['x89']] = float(943.17125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x89']),float(943.17125000)))
        if('y89' in ix_):
            self.x0[ix_['y89']] = float(51011.068905)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y89']),float(51011.068905)))
        if('w90' in ix_):
            self.x0[ix_['w90']] = float(7.4100000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w90']),float(7.4100000000)))
        if('sigt90' in ix_):
            self.x0[ix_['sigt90']] = float(50.755029735)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt90']),float(50.755029735)))
        if('sigr90' in ix_):
            self.x0[ix_['sigr90']] = float(57.752273609)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr90']),float(57.752273609)))
        if('x90' in ix_):
            self.x0[ix_['x90']] = float(946.88500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x90']),float(946.88500000)))
        if('y90' in ix_):
            self.x0[ix_['y90']] = float(51199.747597)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y90']),float(51199.747597)))
        if('w91' in ix_):
            self.x0[ix_['w91']] = float(7.3750000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w91']),float(7.3750000000)))
        if('sigt91' in ix_):
            self.x0[ix_['sigt91']] = float(50.652260905)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt91']),float(50.652260905)))
        if('sigr91' in ix_):
            self.x0[ix_['sigr91']] = float(57.797128035)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr91']),float(57.797128035)))
        if('x91' in ix_):
            self.x0[ix_['x91']] = float(950.58125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x91']),float(950.58125000)))
        if('y91' in ix_):
            self.x0[ix_['y91']] = float(51387.161395)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y91']),float(51387.161395)))
        if('w92' in ix_):
            self.x0[ix_['w92']] = float(7.3400000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w92']),float(7.3400000000)))
        if('sigt92' in ix_):
            self.x0[ix_['sigt92']] = float(50.547275164)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt92']),float(50.547275164)))
        if('sigr92' in ix_):
            self.x0[ix_['sigr92']] = float(57.841788143)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr92']),float(57.841788143)))
        if('x92' in ix_):
            self.x0[ix_['x92']] = float(954.26000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x92']),float(954.26000000)))
        if('y92' in ix_):
            self.x0[ix_['y92']] = float(51573.305751)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y92']),float(51573.305751)))
        if('w93' in ix_):
            self.x0[ix_['w93']] = float(7.3050000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w93']),float(7.3050000000)))
        if('sigt93' in ix_):
            self.x0[ix_['sigt93']] = float(50.440046782)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt93']),float(50.440046782)))
        if('sigr93' in ix_):
            self.x0[ix_['sigr93']] = float(57.886262921)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr93']),float(57.886262921)))
        if('x93' in ix_):
            self.x0[ix_['x93']] = float(957.92125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x93']),float(957.92125000)))
        if('y93' in ix_):
            self.x0[ix_['y93']] = float(51758.176137)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y93']),float(51758.176137)))
        if('w94' in ix_):
            self.x0[ix_['w94']] = float(7.2700000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w94']),float(7.2700000000)))
        if('sigt94' in ix_):
            self.x0[ix_['sigt94']] = float(50.330550025)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt94']),float(50.330550025)))
        if('sigr94' in ix_):
            self.x0[ix_['sigr94']] = float(57.930561476)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr94']),float(57.930561476)))
        if('x94' in ix_):
            self.x0[ix_['x94']] = float(961.56500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x94']),float(961.56500000)))
        if('y94' in ix_):
            self.x0[ix_['y94']] = float(51941.768047)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y94']),float(51941.768047)))
        if('w95' in ix_):
            self.x0[ix_['w95']] = float(7.2350000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w95']),float(7.2350000000)))
        if('sigt95' in ix_):
            self.x0[ix_['sigt95']] = float(50.218759151)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt95']),float(50.218759151)))
        if('sigr95' in ix_):
            self.x0[ix_['sigr95']] = float(57.974693041)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr95']),float(57.974693041)))
        if('x95' in ix_):
            self.x0[ix_['x95']] = float(965.19125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x95']),float(965.19125000)))
        if('y95' in ix_):
            self.x0[ix_['y95']] = float(52124.077002)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y95']),float(52124.077002)))
        if('w96' in ix_):
            self.x0[ix_['w96']] = float(7.2000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w96']),float(7.2000000000)))
        if('sigt96' in ix_):
            self.x0[ix_['sigt96']] = float(50.104648409)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt96']),float(50.104648409)))
        if('sigr96' in ix_):
            self.x0[ix_['sigr96']] = float(58.018666983)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr96']),float(58.018666983)))
        if('x96' in ix_):
            self.x0[ix_['x96']] = float(968.80000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x96']),float(968.80000000)))
        if('y96' in ix_):
            self.x0[ix_['y96']] = float(52305.098550)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y96']),float(52305.098550)))
        if('w97' in ix_):
            self.x0[ix_['w97']] = float(7.1712500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w97']),float(7.1712500000)))
        if('sigt97' in ix_):
            self.x0[ix_['sigt97']] = float(49.972881021)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt97']),float(49.972881021)))
        if('sigr97' in ix_):
            self.x0[ix_['sigr97']] = float(58.011883331)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr97']),float(58.011883331)))
        if('x97' in ix_):
            self.x0[ix_['x97']] = float(972.39281250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x97']),float(972.39281250)))
        if('y97' in ix_):
            self.x0[ix_['y97']] = float(52484.878923)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y97']),float(52484.878923)))
        if('w98' in ix_):
            self.x0[ix_['w98']] = float(7.1425000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w98']),float(7.1425000000)))
        if('sigt98' in ix_):
            self.x0[ix_['sigt98']] = float(49.838337287)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt98']),float(49.838337287)))
        if('sigr98' in ix_):
            self.x0[ix_['sigr98']] = float(58.004462269)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr98']),float(58.004462269)))
        if('x98' in ix_):
            self.x0[ix_['x98']] = float(975.97125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x98']),float(975.97125000)))
        if('y98' in ix_):
            self.x0[ix_['y98']] = float(52663.463510)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y98']),float(52663.463510)))
        if('w99' in ix_):
            self.x0[ix_['w99']] = float(7.1137500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w99']),float(7.1137500000)))
        if('sigt99' in ix_):
            self.x0[ix_['sigt99']] = float(49.700989951)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt99']),float(49.700989951)))
        if('sigr99' in ix_):
            self.x0[ix_['sigr99']] = float(57.996401618)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr99']),float(57.996401618)))
        if('x99' in ix_):
            self.x0[ix_['x99']] = float(979.53531250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x99']),float(979.53531250)))
        if('y99' in ix_):
            self.x0[ix_['y99']] = float(52840.846195)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y99']),float(52840.846195)))
        if('w100' in ix_):
            self.x0[ix_['w100']] = float(7.0850000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w100']),float(7.0850000000)))
        if('sigt100' in ix_):
            self.x0[ix_['sigt100']] = float(49.560811600)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt100']),float(49.560811600)))
        if('sigr100' in ix_):
            self.x0[ix_['sigr100']] = float(57.987699220)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr100']),float(57.987699220)))
        if('x100' in ix_):
            self.x0[ix_['x100']] = float(983.08500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x100']),float(983.08500000)))
        if('y100' in ix_):
            self.x0[ix_['y100']] = float(53017.020887)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y100']),float(53017.020887)))
        if('w101' in ix_):
            self.x0[ix_['w101']] = float(7.0562500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w101']),float(7.0562500000)))
        if('sigt101' in ix_):
            self.x0[ix_['sigt101']] = float(49.417774661)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt101']),float(49.417774661)))
        if('sigr101' in ix_):
            self.x0[ix_['sigr101']] = float(57.978352944)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr101']),float(57.978352944)))
        if('x101' in ix_):
            self.x0[ix_['x101']] = float(986.62031250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x101']),float(986.62031250)))
        if('y101' in ix_):
            self.x0[ix_['y101']] = float(53191.981517)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y101']),float(53191.981517)))
        if('w102' in ix_):
            self.x0[ix_['w102']] = float(7.0275000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w102']),float(7.0275000000)))
        if('sigt102' in ix_):
            self.x0[ix_['sigt102']] = float(49.271851404)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt102']),float(49.271851404)))
        if('sigr102' in ix_):
            self.x0[ix_['sigr102']] = float(57.968360681)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr102']),float(57.968360681)))
        if('x102' in ix_):
            self.x0[ix_['x102']] = float(990.14125000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x102']),float(990.14125000)))
        if('y102' in ix_):
            self.x0[ix_['y102']] = float(53365.722044)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y102']),float(53365.722044)))
        if('w103' in ix_):
            self.x0[ix_['w103']] = float(6.9987500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w103']),float(6.9987500000)))
        if('sigt103' in ix_):
            self.x0[ix_['sigt103']] = float(49.123013943)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt103']),float(49.123013943)))
        if('sigr103' in ix_):
            self.x0[ix_['sigr103']] = float(57.957720351)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr103']),float(57.957720351)))
        if('x103' in ix_):
            self.x0[ix_['x103']] = float(993.64781250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x103']),float(993.64781250)))
        if('y103' in ix_):
            self.x0[ix_['y103']] = float(53538.236452)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y103']),float(53538.236452)))
        if('w104' in ix_):
            self.x0[ix_['w104']] = float(6.9700000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w104']),float(6.9700000000)))
        if('sigt104' in ix_):
            self.x0[ix_['sigt104']] = float(48.971234234)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt104']),float(48.971234234)))
        if('sigr104' in ix_):
            self.x0[ix_['sigr104']] = float(57.946429896)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr104']),float(57.946429896)))
        if('x104' in ix_):
            self.x0[ix_['x104']] = float(997.14000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x104']),float(997.14000000)))
        if('y104' in ix_):
            self.x0[ix_['y104']] = float(53709.518751)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y104']),float(53709.518751)))
        if('w105' in ix_):
            self.x0[ix_['w105']] = float(6.9412500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w105']),float(6.9412500000)))
        if('sigt105' in ix_):
            self.x0[ix_['sigt105']] = float(48.816484081)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt105']),float(48.816484081)))
        if('sigr105' in ix_):
            self.x0[ix_['sigr105']] = float(57.934487286)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr105']),float(57.934487286)))
        if('x105' in ix_):
            self.x0[ix_['x105']] = float(1000.6178125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x105']),float(1000.6178125)))
        if('y105' in ix_):
            self.x0[ix_['y105']] = float(53879.562982)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y105']),float(53879.562982)))
        if('w106' in ix_):
            self.x0[ix_['w106']] = float(6.9125000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w106']),float(6.9125000000)))
        if('sigt106' in ix_):
            self.x0[ix_['sigt106']] = float(48.658735129)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt106']),float(48.658735129)))
        if('sigr106' in ix_):
            self.x0[ix_['sigr106']] = float(57.921890519)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr106']),float(57.921890519)))
        if('x106' in ix_):
            self.x0[ix_['x106']] = float(1004.0812500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x106']),float(1004.0812500)))
        if('y106' in ix_):
            self.x0[ix_['y106']] = float(54048.363213)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y106']),float(54048.363213)))
        if('w107' in ix_):
            self.x0[ix_['w107']] = float(6.8837500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w107']),float(6.8837500000)))
        if('sigt107' in ix_):
            self.x0[ix_['sigt107']] = float(48.497958873)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt107']),float(48.497958873)))
        if('sigr107' in ix_):
            self.x0[ix_['sigr107']] = float(57.908637618)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr107']),float(57.908637618)))
        if('x107' in ix_):
            self.x0[ix_['x107']] = float(1007.5303125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x107']),float(1007.5303125)))
        if('y107' in ix_):
            self.x0[ix_['y107']] = float(54215.913546)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y107']),float(54215.913546)))
        if('w108' in ix_):
            self.x0[ix_['w108']] = float(6.8550000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w108']),float(6.8550000000)))
        if('sigt108' in ix_):
            self.x0[ix_['sigt108']] = float(48.334126653)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt108']),float(48.334126653)))
        if('sigr108' in ix_):
            self.x0[ix_['sigr108']] = float(57.894726636)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr108']),float(57.894726636)))
        if('x108' in ix_):
            self.x0[ix_['x108']] = float(1010.9650000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x108']),float(1010.9650000)))
        if('y108' in ix_):
            self.x0[ix_['y108']] = float(54382.208112)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y108']),float(54382.208112)))
        if('w109' in ix_):
            self.x0[ix_['w109']] = float(6.8262500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w109']),float(6.8262500000)))
        if('sigt109' in ix_):
            self.x0[ix_['sigt109']] = float(48.167209655)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt109']),float(48.167209655)))
        if('sigr109' in ix_):
            self.x0[ix_['sigr109']] = float(57.880155655)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr109']),float(57.880155655)))
        if('x109' in ix_):
            self.x0[ix_['x109']] = float(1014.3853125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x109']),float(1014.3853125)))
        if('y109' in ix_):
            self.x0[ix_['y109']] = float(54547.241075)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y109']),float(54547.241075)))
        if('w110' in ix_):
            self.x0[ix_['w110']] = float(6.7975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w110']),float(6.7975000000)))
        if('sigt110' in ix_):
            self.x0[ix_['sigt110']] = float(47.997178915)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt110']),float(47.997178915)))
        if('sigr110' in ix_):
            self.x0[ix_['sigr110']] = float(57.864922786)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr110']),float(57.864922786)))
        if('x110' in ix_):
            self.x0[ix_['x110']] = float(1017.7912500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x110']),float(1017.7912500)))
        if('y110' in ix_):
            self.x0[ix_['y110']] = float(54711.006635)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y110']),float(54711.006635)))
        if('w111' in ix_):
            self.x0[ix_['w111']] = float(6.7687500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w111']),float(6.7687500000)))
        if('sigt111' in ix_):
            self.x0[ix_['sigt111']] = float(47.824005317)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt111']),float(47.824005317)))
        if('sigr111' in ix_):
            self.x0[ix_['sigr111']] = float(57.849026171)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr111']),float(57.849026171)))
        if('x111' in ix_):
            self.x0[ix_['x111']] = float(1021.1828125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x111']),float(1021.1828125)))
        if('y111' in ix_):
            self.x0[ix_['y111']] = float(54873.499025)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y111']),float(54873.499025)))
        if('w112' in ix_):
            self.x0[ix_['w112']] = float(6.7400000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w112']),float(6.7400000000)))
        if('sigt112' in ix_):
            self.x0[ix_['sigt112']] = float(47.647659596)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt112']),float(47.647659596)))
        if('sigr112' in ix_):
            self.x0[ix_['sigr112']] = float(57.832463983)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr112']),float(57.832463983)))
        if('x112' in ix_):
            self.x0[ix_['x112']] = float(1024.5600000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x112']),float(1024.5600000)))
        if('y112' in ix_):
            self.x0[ix_['y112']] = float(55034.712515)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y112']),float(55034.712515)))
        if('w113' in ix_):
            self.x0[ix_['w113']] = float(6.7187500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w113']),float(6.7187500000)))
        if('sigt113' in ix_):
            self.x0[ix_['sigt113']] = float(47.448591006)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt113']),float(47.448591006)))
        if('sigr113' in ix_):
            self.x0[ix_['sigr113']] = float(57.750663878)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr113']),float(57.750663878)))
        if('x113' in ix_):
            self.x0[ix_['x113']] = float(1027.9246875)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x113']),float(1027.9246875)))
        if('y113' in ix_):
            self.x0[ix_['y113']] = float(55194.697627)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y113']),float(55194.697627)))
        if('w114' in ix_):
            self.x0[ix_['w114']] = float(6.6975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w114']),float(6.6975000000)))
        if('sigt114' in ix_):
            self.x0[ix_['sigt114']] = float(47.245860884)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt114']),float(47.245860884)))
        if('sigr114' in ix_):
            self.x0[ix_['sigr114']] = float(57.667755969)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr114']),float(57.667755969)))
        if('x114' in ix_):
            self.x0[ix_['x114']] = float(1031.2787500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x114']),float(1031.2787500)))
        if('y114' in ix_):
            self.x0[ix_['y114']] = float(55353.503720)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y114']),float(55353.503720)))
        if('w115' in ix_):
            self.x0[ix_['w115']] = float(6.6762500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w115']),float(6.6762500000)))
        if('sigt115' in ix_):
            self.x0[ix_['sigt115']] = float(47.039439651)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt115']),float(47.039439651)))
        if('sigr115' in ix_):
            self.x0[ix_['sigr115']] = float(57.583729120)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr115']),float(57.583729120)))
        if('x115' in ix_):
            self.x0[ix_['x115']] = float(1034.6221875)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x115']),float(1034.6221875)))
        if('y115' in ix_):
            self.x0[ix_['y115']] = float(55511.122773)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y115']),float(55511.122773)))
        if('w116' in ix_):
            self.x0[ix_['w116']] = float(6.6550000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w116']),float(6.6550000000)))
        if('sigt116' in ix_):
            self.x0[ix_['sigt116']] = float(46.829297460)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt116']),float(46.829297460)))
        if('sigr116' in ix_):
            self.x0[ix_['sigr116']] = float(57.498572198)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr116']),float(57.498572198)))
        if('x116' in ix_):
            self.x0[ix_['x116']] = float(1037.9550000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x116']),float(1037.9550000)))
        if('y116' in ix_):
            self.x0[ix_['y116']] = float(55667.546782)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y116']),float(55667.546782)))
        if('w117' in ix_):
            self.x0[ix_['w117']] = float(6.6337500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w117']),float(6.6337500000)))
        if('sigt117' in ix_):
            self.x0[ix_['sigt117']] = float(46.615404203)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt117']),float(46.615404203)))
        if('sigr117' in ix_):
            self.x0[ix_['sigr117']] = float(57.412274068)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr117']),float(57.412274068)))
        if('x117' in ix_):
            self.x0[ix_['x117']] = float(1041.2771875)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x117']),float(1041.2771875)))
        if('y117' in ix_):
            self.x0[ix_['y117']] = float(55822.767760)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y117']),float(55822.767760)))
        if('w118' in ix_):
            self.x0[ix_['w118']] = float(6.6125000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w118']),float(6.6125000000)))
        if('sigt118' in ix_):
            self.x0[ix_['sigt118']] = float(46.397729510)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt118']),float(46.397729510)))
        if('sigr118' in ix_):
            self.x0[ix_['sigr118']] = float(57.324823593)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr118']),float(57.324823593)))
        if('x118' in ix_):
            self.x0[ix_['x118']] = float(1044.5887500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x118']),float(1044.5887500)))
        if('y118' in ix_):
            self.x0[ix_['y118']] = float(55976.777741)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y118']),float(55976.777741)))
        if('w119' in ix_):
            self.x0[ix_['w119']] = float(6.5912500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w119']),float(6.5912500000)))
        if('sigt119' in ix_):
            self.x0[ix_['sigt119']] = float(46.176242754)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt119']),float(46.176242754)))
        if('sigr119' in ix_):
            self.x0[ix_['sigr119']] = float(57.236209626)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr119']),float(57.236209626)))
        if('x119' in ix_):
            self.x0[ix_['x119']] = float(1047.8896875)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x119']),float(1047.8896875)))
        if('y119' in ix_):
            self.x0[ix_['y119']] = float(56129.568777)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y119']),float(56129.568777)))
        if('w120' in ix_):
            self.x0[ix_['w120']] = float(6.5700000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w120']),float(6.5700000000)))
        if('sigt120' in ix_):
            self.x0[ix_['sigt120']] = float(45.950913048)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt120']),float(45.950913048)))
        if('sigr120' in ix_):
            self.x0[ix_['sigr120']] = float(57.146421013)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr120']),float(57.146421013)))
        if('x120' in ix_):
            self.x0[ix_['x120']] = float(1051.1800000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x120']),float(1051.1800000)))
        if('y120' in ix_):
            self.x0[ix_['y120']] = float(56281.132942)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y120']),float(56281.132942)))
        if('w121' in ix_):
            self.x0[ix_['w121']] = float(6.5412500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w121']),float(6.5412500000)))
        if('sigt121' in ix_):
            self.x0[ix_['sigt121']] = float(45.741494802)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt121']),float(45.741494802)))
        if('sigr121' in ix_):
            self.x0[ix_['sigr121']] = float(57.120910879)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr121']),float(57.120910879)))
        if('x121' in ix_):
            self.x0[ix_['x121']] = float(1054.4578125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x121']),float(1054.4578125)))
        if('y121' in ix_):
            self.x0[ix_['y121']] = float(56431.408955)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y121']),float(56431.408955)))
        if('w122' in ix_):
            self.x0[ix_['w122']] = float(6.5125000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w122']),float(6.5125000000)))
        if('sigt122' in ix_):
            self.x0[ix_['sigt122']] = float(45.528600901)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt122']),float(45.528600901)))
        if('sigr122' in ix_):
            self.x0[ix_['sigr122']] = float(57.094665862)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr122']),float(57.094665862)))
        if('x122' in ix_):
            self.x0[ix_['x122']] = float(1057.7212500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x122']),float(1057.7212500)))
        if('y122' in ix_):
            self.x0[ix_['y122']] = float(56580.336846)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y122']),float(56580.336846)))
        if('w123' in ix_):
            self.x0[ix_['w123']] = float(6.4837500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w123']),float(6.4837500000)))
        if('sigt123' in ix_):
            self.x0[ix_['sigt123']] = float(45.312199922)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt123']),float(45.312199922)))
        if('sigr123' in ix_):
            self.x0[ix_['sigr123']] = float(57.067684218)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr123']),float(57.067684218)))
        if('x123' in ix_):
            self.x0[ix_['x123']] = float(1060.9703125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x123']),float(1060.9703125)))
        if('y123' in ix_):
            self.x0[ix_['y123']] = float(56727.911344)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y123']),float(56727.911344)))
        if('w124' in ix_):
            self.x0[ix_['w124']] = float(6.4550000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w124']),float(6.4550000000)))
        if('sigt124' in ix_):
            self.x0[ix_['sigt124']] = float(45.092260303)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt124']),float(45.092260303)))
        if('sigr124' in ix_):
            self.x0[ix_['sigr124']] = float(57.039964224)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr124']),float(57.039964224)))
        if('x124' in ix_):
            self.x0[ix_['x124']] = float(1064.2050000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x124']),float(1064.2050000)))
        if('y124' in ix_):
            self.x0[ix_['y124']] = float(56874.127223)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y124']),float(56874.127223)))
        if('w125' in ix_):
            self.x0[ix_['w125']] = float(6.4262500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w125']),float(6.4262500000)))
        if('sigt125' in ix_):
            self.x0[ix_['sigt125']] = float(44.868750344)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt125']),float(44.868750344)))
        if('sigr125' in ix_):
            self.x0[ix_['sigr125']] = float(57.011504190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr125']),float(57.011504190)))
        if('x125' in ix_):
            self.x0[ix_['x125']] = float(1067.4253125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x125']),float(1067.4253125)))
        if('y125' in ix_):
            self.x0[ix_['y125']] = float(57018.979310)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y125']),float(57018.979310)))
        if('w126' in ix_):
            self.x0[ix_['w126']] = float(6.3975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w126']),float(6.3975000000)))
        if('sigt126' in ix_):
            self.x0[ix_['sigt126']] = float(44.641638205)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt126']),float(44.641638205)))
        if('sigr126' in ix_):
            self.x0[ix_['sigr126']] = float(56.982302451)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr126']),float(56.982302451)))
        if('x126' in ix_):
            self.x0[ix_['x126']] = float(1070.6312500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x126']),float(1070.6312500)))
        if('y126' in ix_):
            self.x0[ix_['y126']] = float(57162.462482)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y126']),float(57162.462482)))
        if('w127' in ix_):
            self.x0[ix_['w127']] = float(6.3687500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w127']),float(6.3687500000)))
        if('sigt127' in ix_):
            self.x0[ix_['sigt127']] = float(44.410891909)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt127']),float(44.410891909)))
        if('sigr127' in ix_):
            self.x0[ix_['sigr127']] = float(56.952357376)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr127']),float(56.952357376)))
        if('x127' in ix_):
            self.x0[ix_['x127']] = float(1073.8228125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x127']),float(1073.8228125)))
        if('y127' in ix_):
            self.x0[ix_['y127']] = float(57304.571669)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y127']),float(57304.571669)))
        if('w128' in ix_):
            self.x0[ix_['w128']] = float(6.3400000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w128']),float(6.3400000000)))
        if('sigt128' in ix_):
            self.x0[ix_['sigt128']] = float(44.176479340)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt128']),float(44.176479340)))
        if('sigr128' in ix_):
            self.x0[ix_['sigr128']] = float(56.921667364)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr128']),float(56.921667364)))
        if('x128' in ix_):
            self.x0[ix_['x128']] = float(1077.0000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x128']),float(1077.0000000)))
        if('y128' in ix_):
            self.x0[ix_['y128']] = float(57445.301855)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y128']),float(57445.301855)))
        if('w129' in ix_):
            self.x0[ix_['w129']] = float(6.3162500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w129']),float(6.3162500000)))
        if('sigt129' in ix_):
            self.x0[ix_['sigt129']] = float(43.924748683)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt129']),float(43.924748683)))
        if('sigr129' in ix_):
            self.x0[ix_['sigr129']] = float(56.845155310)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr129']),float(56.845155310)))
        if('x129' in ix_):
            self.x0[ix_['x129']] = float(1080.1640625)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x129']),float(1080.1640625)))
        if('y129' in ix_):
            self.x0[ix_['y129']] = float(57584.681499)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y129']),float(57584.681499)))
        if('w130' in ix_):
            self.x0[ix_['w130']] = float(6.2925000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w130']),float(6.2925000000)))
        if('sigt130' in ix_):
            self.x0[ix_['sigt130']] = float(43.668982053)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt130']),float(43.668982053)))
        if('sigr130' in ix_):
            self.x0[ix_['sigr130']] = float(56.767522016)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr130']),float(56.767522016)))
        if('x130' in ix_):
            self.x0[ix_['x130']] = float(1083.3162500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x130']),float(1083.3162500)))
        if('y130' in ix_):
            self.x0[ix_['y130']] = float(57722.738189)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y130']),float(57722.738189)))
        if('w131' in ix_):
            self.x0[ix_['w131']] = float(6.2687500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w131']),float(6.2687500000)))
        if('sigt131' in ix_):
            self.x0[ix_['sigt131']] = float(43.409146055)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt131']),float(43.409146055)))
        if('sigr131' in ix_):
            self.x0[ix_['sigr131']] = float(56.688758490)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr131']),float(56.688758490)))
        if('x131' in ix_):
            self.x0[ix_['x131']] = float(1086.4565625)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x131']),float(1086.4565625)))
        if('y131' in ix_):
            self.x0[ix_['y131']] = float(57859.465228)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y131']),float(57859.465228)))
        if('w132' in ix_):
            self.x0[ix_['w132']] = float(6.2450000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w132']),float(6.2450000000)))
        if('sigt132' in ix_):
            self.x0[ix_['sigt132']] = float(43.145207079)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt132']),float(43.145207079)))
        if('sigr132' in ix_):
            self.x0[ix_['sigr132']] = float(56.608855703)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr132']),float(56.608855703)))
        if('x132' in ix_):
            self.x0[ix_['x132']] = float(1089.5850000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x132']),float(1089.5850000)))
        if('y132' in ix_):
            self.x0[ix_['y132']] = float(57994.855954)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y132']),float(57994.855954)))
        if('w133' in ix_):
            self.x0[ix_['w133']] = float(6.2212500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w133']),float(6.2212500000)))
        if('sigt133' in ix_):
            self.x0[ix_['sigt133']] = float(42.877131300)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt133']),float(42.877131300)))
        if('sigr133' in ix_):
            self.x0[ix_['sigr133']] = float(56.527804595)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr133']),float(56.527804595)))
        if('x133' in ix_):
            self.x0[ix_['x133']] = float(1092.7015625)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x133']),float(1092.7015625)))
        if('y133' in ix_):
            self.x0[ix_['y133']] = float(58128.903746)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y133']),float(58128.903746)))
        if('w134' in ix_):
            self.x0[ix_['w134']] = float(6.1975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w134']),float(6.1975000000)))
        if('sigt134' in ix_):
            self.x0[ix_['sigt134']] = float(42.604884678)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt134']),float(42.604884678)))
        if('sigr134' in ix_):
            self.x0[ix_['sigr134']] = float(56.445596072)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr134']),float(56.445596072)))
        if('x134' in ix_):
            self.x0[ix_['x134']] = float(1095.8062500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x134']),float(1095.8062500)))
        if('y134' in ix_):
            self.x0[ix_['y134']] = float(58261.602028)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y134']),float(58261.602028)))
        if('w135' in ix_):
            self.x0[ix_['w135']] = float(6.1737500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w135']),float(6.1737500000)))
        if('sigt135' in ix_):
            self.x0[ix_['sigt135']] = float(42.328432957)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt135']),float(42.328432957)))
        if('sigr135' in ix_):
            self.x0[ix_['sigr135']] = float(56.362220999)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr135']),float(56.362220999)))
        if('x135' in ix_):
            self.x0[ix_['x135']] = float(1098.8990625)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x135']),float(1098.8990625)))
        if('y135' in ix_):
            self.x0[ix_['y135']] = float(58392.944262)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y135']),float(58392.944262)))
        if('w136' in ix_):
            self.x0[ix_['w136']] = float(6.1500000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w136']),float(6.1500000000)))
        if('sigt136' in ix_):
            self.x0[ix_['sigt136']] = float(42.047741672)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt136']),float(42.047741672)))
        if('sigr136' in ix_):
            self.x0[ix_['sigr136']] = float(56.277670207)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr136']),float(56.277670207)))
        if('x136' in ix_):
            self.x0[ix_['x136']] = float(1101.9800000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x136']),float(1101.9800000)))
        if('y136' in ix_):
            self.x0[ix_['y136']] = float(58522.923955)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y136']),float(58522.923955)))
        if('w137' in ix_):
            self.x0[ix_['w137']] = float(6.1150000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w137']),float(6.1150000000)))
        if('sigt137' in ix_):
            self.x0[ix_['sigt137']] = float(41.794038608)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt137']),float(41.794038608)))
        if('sigr137' in ix_):
            self.x0[ix_['sigr137']] = float(56.295428092)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr137']),float(56.295428092)))
        if('x137' in ix_):
            self.x0[ix_['x137']] = float(1105.0462500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x137']),float(1105.0462500)))
        if('y137' in ix_):
            self.x0[ix_['y137']] = float(58651.464995)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y137']),float(58651.464995)))
        if('w138' in ix_):
            self.x0[ix_['w138']] = float(6.0800000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w138']),float(6.0800000000)))
        if('sigt138' in ix_):
            self.x0[ix_['sigt138']] = float(41.536786267)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt138']),float(41.536786267)))
        if('sigr138' in ix_):
            self.x0[ix_['sigr138']] = float(56.313098515)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr138']),float(56.313098515)))
        if('x138' in ix_):
            self.x0[ix_['x138']] = float(1108.0950000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x138']),float(1108.0950000)))
        if('y138' in ix_):
            self.x0[ix_['y138']] = float(58778.493546)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y138']),float(58778.493546)))
        if('w139' in ix_):
            self.x0[ix_['w139']] = float(6.0450000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w139']),float(6.0450000000)))
        if('sigt139' in ix_):
            self.x0[ix_['sigt139']] = float(41.275954987)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt139']),float(41.275954987)))
        if('sigr139' in ix_):
            self.x0[ix_['sigr139']] = float(56.330696252)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr139']),float(56.330696252)))
        if('x139' in ix_):
            self.x0[ix_['x139']] = float(1111.1262500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x139']),float(1111.1262500)))
        if('y139' in ix_):
            self.x0[ix_['y139']] = float(58904.007748)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y139']),float(58904.007748)))
        if('w140' in ix_):
            self.x0[ix_['w140']] = float(6.0100000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w140']),float(6.0100000000)))
        if('sigt140' in ix_):
            self.x0[ix_['sigt140']] = float(41.011515151)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt140']),float(41.011515151)))
        if('sigr140' in ix_):
            self.x0[ix_['sigr140']] = float(56.348236444)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr140']),float(56.348236444)))
        if('x140' in ix_):
            self.x0[ix_['x140']] = float(1114.1400000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x140']),float(1114.1400000)))
        if('y140' in ix_):
            self.x0[ix_['y140']] = float(59028.005837)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y140']),float(59028.005837)))
        if('w141' in ix_):
            self.x0[ix_['w141']] = float(5.9750000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w141']),float(5.9750000000)))
        if('sigt141' in ix_):
            self.x0[ix_['sigt141']] = float(40.743437190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt141']),float(40.743437190)))
        if('sigr141' in ix_):
            self.x0[ix_['sigr141']] = float(56.365734613)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr141']),float(56.365734613)))
        if('x141' in ix_):
            self.x0[ix_['x141']] = float(1117.1362500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x141']),float(1117.1362500)))
        if('y141' in ix_):
            self.x0[ix_['y141']] = float(59150.486148)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y141']),float(59150.486148)))
        if('w142' in ix_):
            self.x0[ix_['w142']] = float(5.9400000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w142']),float(5.9400000000)))
        if('sigt142' in ix_):
            self.x0[ix_['sigt142']] = float(40.471691591)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt142']),float(40.471691591)))
        if('sigr142' in ix_):
            self.x0[ix_['sigr142']] = float(56.383206675)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr142']),float(56.383206675)))
        if('x142' in ix_):
            self.x0[ix_['x142']] = float(1120.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x142']),float(1120.1150000)))
        if('y142' in ix_):
            self.x0[ix_['y142']] = float(59271.447119)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y142']),float(59271.447119)))
        if('w143' in ix_):
            self.x0[ix_['w143']] = float(5.9050000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w143']),float(5.9050000000)))
        if('sigt143' in ix_):
            self.x0[ix_['sigt143']] = float(40.196248899)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt143']),float(40.196248899)))
        if('sigr143' in ix_):
            self.x0[ix_['sigr143']] = float(56.400668955)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr143']),float(56.400668955)))
        if('x143' in ix_):
            self.x0[ix_['x143']] = float(1123.0762500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x143']),float(1123.0762500)))
        if('y143' in ix_):
            self.x0[ix_['y143']] = float(59390.887293)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y143']),float(59390.887293)))
        if('w144' in ix_):
            self.x0[ix_['w144']] = float(5.8700000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w144']),float(5.8700000000)))
        if('sigt144' in ix_):
            self.x0[ix_['sigt144']] = float(39.917079719)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt144']),float(39.917079719)))
        if('sigr144' in ix_):
            self.x0[ix_['sigr144']] = float(56.418138202)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr144']),float(56.418138202)))
        if('x144' in ix_):
            self.x0[ix_['x144']] = float(1126.0200000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x144']),float(1126.0200000)))
        if('y144' in ix_):
            self.x0[ix_['y144']] = float(59508.805320)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y144']),float(59508.805320)))
        if('w145' in ix_):
            self.x0[ix_['w145']] = float(5.8412500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w145']),float(5.8412500000)))
        if('sigt145' in ix_):
            self.x0[ix_['sigt145']] = float(39.615894752)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt145']),float(39.615894752)))
        if('sigr145' in ix_):
            self.x0[ix_['sigr145']] = float(56.375167857)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr145']),float(56.375167857)))
        if('x145' in ix_):
            self.x0[ix_['x145']] = float(1128.9478125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x145']),float(1128.9478125)))
        if('y145' in ix_):
            self.x0[ix_['y145']] = float(59625.235221)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y145']),float(59625.235221)))
        if('w146' in ix_):
            self.x0[ix_['w146']] = float(5.8125000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w146']),float(5.8125000000)))
        if('sigt146' in ix_):
            self.x0[ix_['sigt146']] = float(39.310443929)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt146']),float(39.310443929)))
        if('sigr146' in ix_):
            self.x0[ix_['sigr146']] = float(56.331441829)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr146']),float(56.331441829)))
        if('x146' in ix_):
            self.x0[ix_['x146']] = float(1131.8612500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x146']),float(1131.8612500)))
        if('y146' in ix_):
            self.x0[ix_['y146']] = float(59740.209796)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y146']),float(59740.209796)))
        if('w147' in ix_):
            self.x0[ix_['w147']] = float(5.7837500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w147']),float(5.7837500000)))
        if('sigt147' in ix_):
            self.x0[ix_['sigt147']] = float(39.000692735)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt147']),float(39.000692735)))
        if('sigr147' in ix_):
            self.x0[ix_['sigr147']] = float(56.286959540)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr147']),float(56.286959540)))
        if('x147' in ix_):
            self.x0[ix_['x147']] = float(1134.7603125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x147']),float(1134.7603125)))
        if('y147' in ix_):
            self.x0[ix_['y147']] = float(59853.725349)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y147']),float(59853.725349)))
        if('w148' in ix_):
            self.x0[ix_['w148']] = float(5.7550000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w148']),float(5.7550000000)))
        if('sigt148' in ix_):
            self.x0[ix_['sigt148']] = float(38.686606531)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt148']),float(38.686606531)))
        if('sigr148' in ix_):
            self.x0[ix_['sigr148']] = float(56.241720481)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr148']),float(56.241720481)))
        if('x148' in ix_):
            self.x0[ix_['x148']] = float(1137.6450000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x148']),float(1137.6450000)))
        if('y148' in ix_):
            self.x0[ix_['y148']] = float(59965.778269)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y148']),float(59965.778269)))
        if('w149' in ix_):
            self.x0[ix_['w149']] = float(5.7262500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w149']),float(5.7262500000)))
        if('sigt149' in ix_):
            self.x0[ix_['sigt149']] = float(38.368150555)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt149']),float(38.368150555)))
        if('sigr149' in ix_):
            self.x0[ix_['sigr149']] = float(56.195724220)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr149']),float(56.195724220)))
        if('x149' in ix_):
            self.x0[ix_['x149']] = float(1140.5153125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x149']),float(1140.5153125)))
        if('y149' in ix_):
            self.x0[ix_['y149']] = float(60076.365029)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y149']),float(60076.365029)))
        if('w150' in ix_):
            self.x0[ix_['w150']] = float(5.6975000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w150']),float(5.6975000000)))
        if('sigt150' in ix_):
            self.x0[ix_['sigt150']] = float(38.045289920)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt150']),float(38.045289920)))
        if('sigr150' in ix_):
            self.x0[ix_['sigr150']] = float(56.148970399)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr150']),float(56.148970399)))
        if('x150' in ix_):
            self.x0[ix_['x150']] = float(1143.3712500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x150']),float(1143.3712500)))
        if('y150' in ix_):
            self.x0[ix_['y150']] = float(60185.482195)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y150']),float(60185.482195)))
        if('w151' in ix_):
            self.x0[ix_['w151']] = float(5.6687500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w151']),float(5.6687500000)))
        if('sigt151' in ix_):
            self.x0[ix_['sigt151']] = float(37.717989622)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt151']),float(37.717989622)))
        if('sigr151' in ix_):
            self.x0[ix_['sigr151']] = float(56.101458742)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr151']),float(56.101458742)))
        if('x151' in ix_):
            self.x0[ix_['x151']] = float(1146.2128125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x151']),float(1146.2128125)))
        if('y151' in ix_):
            self.x0[ix_['y151']] = float(60293.126418)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y151']),float(60293.126418)))
        if('w152' in ix_):
            self.x0[ix_['w152']] = float(5.6400000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w152']),float(5.6400000000)))
        if('sigt152' in ix_):
            self.x0[ix_['sigt152']] = float(37.386214534)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt152']),float(37.386214534)))
        if('sigr152' in ix_):
            self.x0[ix_['sigr152']] = float(56.053189053)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr152']),float(56.053189053)))
        if('x152' in ix_):
            self.x0[ix_['x152']] = float(1149.0400000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x152']),float(1149.0400000)))
        if('y152' in ix_):
            self.x0[ix_['y152']] = float(60399.294444)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y152']),float(60399.294444)))
        if('w153' in ix_):
            self.x0[ix_['w153']] = float(5.6112500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w153']),float(5.6112500000)))
        if('sigt153' in ix_):
            self.x0[ix_['sigt153']] = float(37.049929410)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt153']),float(37.049929410)))
        if('sigr153' in ix_):
            self.x0[ix_['sigr153']] = float(56.004161223)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr153']),float(56.004161223)))
        if('x153' in ix_):
            self.x0[ix_['x153']] = float(1151.8528125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x153']),float(1151.8528125)))
        if('y153' in ix_):
            self.x0[ix_['y153']] = float(60503.983110)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y153']),float(60503.983110)))
        if('w154' in ix_):
            self.x0[ix_['w154']] = float(5.5825000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w154']),float(5.5825000000)))
        if('sigt154' in ix_):
            self.x0[ix_['sigt154']] = float(36.709098885)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt154']),float(36.709098885)))
        if('sigr154' in ix_):
            self.x0[ix_['sigr154']] = float(55.954375227)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr154']),float(55.954375227)))
        if('x154' in ix_):
            self.x0[ix_['x154']] = float(1154.6512500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x154']),float(1154.6512500)))
        if('y154' in ix_):
            self.x0[ix_['y154']] = float(60607.189351)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y154']),float(60607.189351)))
        if('w155' in ix_):
            self.x0[ix_['w155']] = float(5.5537500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w155']),float(5.5537500000)))
        if('sigt155' in ix_):
            self.x0[ix_['sigt155']] = float(36.363687480)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt155']),float(36.363687480)))
        if('sigr155' in ix_):
            self.x0[ix_['sigr155']] = float(55.903831135)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr155']),float(55.903831135)))
        if('x155' in ix_):
            self.x0[ix_['x155']] = float(1157.4353125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x155']),float(1157.4353125)))
        if('y155' in ix_):
            self.x0[ix_['y155']] = float(60708.910194)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y155']),float(60708.910194)))
        if('w156' in ix_):
            self.x0[ix_['w156']] = float(5.5250000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w156']),float(5.5250000000)))
        if('sigt156' in ix_):
            self.x0[ix_['sigt156']] = float(36.013659597)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt156']),float(36.013659597)))
        if('sigr156' in ix_):
            self.x0[ix_['sigr156']] = float(55.852529108)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr156']),float(55.852529108)))
        if('x156' in ix_):
            self.x0[ix_['x156']] = float(1160.2050000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x156']),float(1160.2050000)))
        if('y156' in ix_):
            self.x0[ix_['y156']] = float(60809.142769)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y156']),float(60809.142769)))
        if('w157' in ix_):
            self.x0[ix_['w157']] = float(5.4962500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w157']),float(5.4962500000)))
        if('sigt157' in ix_):
            self.x0[ix_['sigt157']] = float(35.658979526)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt157']),float(35.658979526)))
        if('sigr157' in ix_):
            self.x0[ix_['sigr157']] = float(55.800469405)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr157']),float(55.800469405)))
        if('x157' in ix_):
            self.x0[ix_['x157']] = float(1162.9603125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x157']),float(1162.9603125)))
        if('y157' in ix_):
            self.x0[ix_['y157']] = float(60907.884303)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y157']),float(60907.884303)))
        if('w158' in ix_):
            self.x0[ix_['w158']] = float(5.4675000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w158']),float(5.4675000000)))
        if('sigt158' in ix_):
            self.x0[ix_['sigt158']] = float(35.299611440)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt158']),float(35.299611440)))
        if('sigr158' in ix_):
            self.x0[ix_['sigr158']] = float(55.747652383)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr158']),float(55.747652383)))
        if('x158' in ix_):
            self.x0[ix_['x158']] = float(1165.7012500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x158']),float(1165.7012500)))
        if('y158' in ix_):
            self.x0[ix_['y158']] = float(61005.132126)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y158']),float(61005.132126)))
        if('w159' in ix_):
            self.x0[ix_['w159']] = float(5.4387500000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w159']),float(5.4387500000)))
        if('sigt159' in ix_):
            self.x0[ix_['sigt159']] = float(34.935519403)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt159']),float(34.935519403)))
        if('sigr159' in ix_):
            self.x0[ix_['sigr159']] = float(55.694078505)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr159']),float(55.694078505)))
        if('x159' in ix_):
            self.x0[ix_['x159']] = float(1168.4278125)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x159']),float(1168.4278125)))
        if('y159' in ix_):
            self.x0[ix_['y159']] = float(61100.883671)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y159']),float(61100.883671)))
        if('w160' in ix_):
            self.x0[ix_['w160']] = float(5.4100000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w160']),float(5.4100000000)))
        if('sigt160' in ix_):
            self.x0[ix_['sigt160']] = float(34.566667368)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt160']),float(34.566667368)))
        if('sigr160' in ix_):
            self.x0[ix_['sigr160']] = float(55.639748340)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr160']),float(55.639748340)))
        if('x160' in ix_):
            self.x0[ix_['x160']] = float(1171.1400000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x160']),float(1171.1400000)))
        if('y160' in ix_):
            self.x0[ix_['y160']] = float(61195.136478)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y160']),float(61195.136478)))
        if('w161' in ix_):
            self.x0[ix_['w161']] = float(5.6025000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w161']),float(5.6025000000)))
        if('sigt161' in ix_):
            self.x0[ix_['sigt161']] = float(33.529237255)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt161']),float(33.529237255)))
        if('sigr161' in ix_):
            self.x0[ix_['sigr161']] = float(53.385743943)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr161']),float(53.385743943)))
        if('x161' in ix_):
            self.x0[ix_['x161']] = float(1173.8931250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x161']),float(1173.8931250)))
        if('y161' in ix_):
            self.x0[ix_['y161']] = float(61288.849783)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y161']),float(61288.849783)))
        if('w162' in ix_):
            self.x0[ix_['w162']] = float(5.7950000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w162']),float(5.7950000000)))
        if('sigt162' in ix_):
            self.x0[ix_['sigt162']] = float(32.521949841)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt162']),float(32.521949841)))
        if('sigr162' in ix_):
            self.x0[ix_['sigr162']] = float(51.273873576)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr162']),float(51.273873576)))
        if('x162' in ix_):
            self.x0[ix_['x162']] = float(1176.7425000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x162']),float(1176.7425000)))
        if('y162' in ix_):
            self.x0[ix_['y162']] = float(61382.927846)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y162']),float(61382.927846)))
        if('w163' in ix_):
            self.x0[ix_['w163']] = float(5.9875000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w163']),float(5.9875000000)))
        if('sigt163' in ix_):
            self.x0[ix_['sigt163']] = float(31.541203081)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt163']),float(31.541203081)))
        if('sigr163' in ix_):
            self.x0[ix_['sigr163']] = float(49.290216376)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr163']),float(49.290216376)))
        if('x163' in ix_):
            self.x0[ix_['x163']] = float(1179.6881250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x163']),float(1179.6881250)))
        if('y163' in ix_):
            self.x0[ix_['y163']] = float(61477.257259)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y163']),float(61477.257259)))
        if('w164' in ix_):
            self.x0[ix_['w164']] = float(6.1800000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w164']),float(6.1800000000)))
        if('sigt164' in ix_):
            self.x0[ix_['sigt164']] = float(30.583856128)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt164']),float(30.583856128)))
        if('sigr164' in ix_):
            self.x0[ix_['sigr164']] = float(47.422585757)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr164']),float(47.422585757)))
        if('x164' in ix_):
            self.x0[ix_['x164']] = float(1182.7300000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x164']),float(1182.7300000)))
        if('y164' in ix_):
            self.x0[ix_['y164']] = float(61571.722555)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y164']),float(61571.722555)))
        if('w165' in ix_):
            self.x0[ix_['w165']] = float(6.8125000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w165']),float(6.8125000000)))
        if('sigt165' in ix_):
            self.x0[ix_['sigt165']] = float(28.755024885)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt165']),float(28.755024885)))
        if('sigr165' in ix_):
            self.x0[ix_['sigr165']] = float(42.704591925)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr165']),float(42.704591925)))
        if('x165' in ix_):
            self.x0[ix_['x165']] = float(1185.9781250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x165']),float(1185.9781250)))
        if('y165' in ix_):
            self.x0[ix_['y165']] = float(61667.948015)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y165']),float(61667.948015)))
        if('w166' in ix_):
            self.x0[ix_['w166']] = float(7.4450000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w166']),float(7.4450000000)))
        if('sigt166' in ix_):
            self.x0[ix_['sigt166']] = float(27.140969642)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt166']),float(27.140969642)))
        if('sigr166' in ix_):
            self.x0[ix_['sigr166']] = float(38.769382165)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr166']),float(38.769382165)))
        if('x166' in ix_):
            self.x0[ix_['x166']] = float(1189.5425000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x166']),float(1189.5425000)))
        if('y166' in ix_):
            self.x0[ix_['y166']] = float(61767.437546)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y166']),float(61767.437546)))
        if('w167' in ix_):
            self.x0[ix_['w167']] = float(8.0775000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w167']),float(8.0775000000)))
        if('sigt167' in ix_):
            self.x0[ix_['sigt167']] = float(25.689067179)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt167']),float(25.689067179)))
        if('sigr167' in ix_):
            self.x0[ix_['sigr167']] = float(35.432574692)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr167']),float(35.432574692)))
        if('x167' in ix_):
            self.x0[ix_['x167']] = float(1193.4231250)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x167']),float(1193.4231250)))
        if('y167' in ix_):
            self.x0[ix_['y167']] = float(61869.829536)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y167']),float(61869.829536)))
        if('w168' in ix_):
            self.x0[ix_['w168']] = float(8.7100000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w168']),float(8.7100000000)))
        if('sigt168' in ix_):
            self.x0[ix_['sigt168']] = float(24.362140259)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt168']),float(24.362140259)))
        if('sigr168' in ix_):
            self.x0[ix_['sigr168']] = float(32.563342979)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr168']),float(32.563342979)))
        if('x168' in ix_):
            self.x0[ix_['x168']] = float(1197.6200000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x168']),float(1197.6200000)))
        if('y168' in ix_):
            self.x0[ix_['y168']] = float(61974.753956)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y168']),float(61974.753956)))
        if('w169' in ix_):
            self.x0[ix_['w169']] = float(9.6750000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w169']),float(9.6750000000)))
        if('sigt169' in ix_):
            self.x0[ix_['sigt169']] = float(22.820212971)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt169']),float(22.820212971)))
        if('sigr169' in ix_):
            self.x0[ix_['sigr169']] = float(29.029268985)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr169']),float(29.029268985)))
        if('x169' in ix_):
            self.x0[ix_['x169']] = float(1202.2162500)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x169']),float(1202.2162500)))
        if('y169' in ix_):
            self.x0[ix_['y169']] = float(62082.998907)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y169']),float(62082.998907)))
        if('w170' in ix_):
            self.x0[ix_['w170']] = float(10.640000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w170']),float(10.640000000)))
        if('sigt170' in ix_):
            self.x0[ix_['sigt170']] = float(21.448594019)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt170']),float(21.448594019)))
        if('sigr170' in ix_):
            self.x0[ix_['sigr170']] = float(26.114653928)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr170']),float(26.114653928)))
        if('x170' in ix_):
            self.x0[ix_['x170']] = float(1207.2950000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x170']),float(1207.2950000)))
        if('y170' in ix_):
            self.x0[ix_['y170']] = float(62195.248557)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y170']),float(62195.248557)))
        if('w171' in ix_):
            self.x0[ix_['w171']] = float(11.820000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w171']),float(11.820000000)))
        if('sigt171' in ix_):
            self.x0[ix_['sigt171']] = float(20.072296567)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt171']),float(20.072296567)))
        if('sigr171' in ix_):
            self.x0[ix_['sigr171']] = float(23.231954054)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr171']),float(23.231954054)))
        if('x171' in ix_):
            self.x0[ix_['x171']] = float(1212.9100000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x171']),float(1212.9100000)))
        if('y171' in ix_):
            self.x0[ix_['y171']] = float(62311.615454)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y171']),float(62311.615454)))
        if('w172' in ix_):
            self.x0[ix_['w172']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w172']),float(13.000000000)))
        if('sigt172' in ix_):
            self.x0[ix_['sigt172']] = float(18.833069241)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt172']),float(18.833069241)))
        if('sigr172' in ix_):
            self.x0[ix_['sigr172']] = float(20.850204821)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr172']),float(20.850204821)))
        if('x172' in ix_):
            self.x0[ix_['x172']] = float(1219.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x172']),float(1219.1150000)))
        if('y172' in ix_):
            self.x0[ix_['y172']] = float(62432.136565)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y172']),float(62432.136565)))
        if('w173' in ix_):
            self.x0[ix_['w173']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w173']),float(13.000000000)))
        if('sigt173' in ix_):
            self.x0[ix_['sigt173']] = float(18.214200054)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt173']),float(18.214200054)))
        if('sigr173' in ix_):
            self.x0[ix_['sigr173']] = float(20.564709344)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr173']),float(20.564709344)))
        if('x173' in ix_):
            self.x0[ix_['x173']] = float(1225.6150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x173']),float(1225.6150000)))
        if('y173' in ix_):
            self.x0[ix_['y173']] = float(62552.540190)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y173']),float(62552.540190)))
        if('w174' in ix_):
            self.x0[ix_['w174']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w174']),float(13.000000000)))
        if('sigt174' in ix_):
            self.x0[ix_['sigt174']] = float(17.589839561)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt174']),float(17.589839561)))
        if('sigr174' in ix_):
            self.x0[ix_['sigr174']] = float(20.276848479)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr174']),float(20.276848479)))
        if('x174' in ix_):
            self.x0[ix_['x174']] = float(1232.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x174']),float(1232.1150000)))
        if('y174' in ix_):
            self.x0[ix_['y174']] = float(62668.903319)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y174']),float(62668.903319)))
        if('w175' in ix_):
            self.x0[ix_['w175']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w175']),float(13.000000000)))
        if('sigt175' in ix_):
            self.x0[ix_['sigt175']] = float(16.959939451)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt175']),float(16.959939451)))
        if('sigr175' in ix_):
            self.x0[ix_['sigr175']] = float(19.986619910)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr175']),float(19.986619910)))
        if('x175' in ix_):
            self.x0[ix_['x175']] = float(1238.6150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x175']),float(1238.6150000)))
        if('y175' in ix_):
            self.x0[ix_['y175']] = float(62781.190101)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y175']),float(62781.190101)))
        if('w176' in ix_):
            self.x0[ix_['w176']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w176']),float(13.000000000)))
        if('sigt176' in ix_):
            self.x0[ix_['sigt176']] = float(16.324451370)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt176']),float(16.324451370)))
        if('sigr176' in ix_):
            self.x0[ix_['sigr176']] = float(19.694021169)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr176']),float(19.694021169)))
        if('x176' in ix_):
            self.x0[ix_['x176']] = float(1245.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x176']),float(1245.1150000)))
        if('y176' in ix_):
            self.x0[ix_['y176']] = float(62889.364371)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y176']),float(62889.364371)))
        if('w177' in ix_):
            self.x0[ix_['w177']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w177']),float(13.000000000)))
        if('sigt177' in ix_):
            self.x0[ix_['sigt177']] = float(15.683326910)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt177']),float(15.683326910)))
        if('sigr177' in ix_):
            self.x0[ix_['sigr177']] = float(19.399049638)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr177']),float(19.399049638)))
        if('x177' in ix_):
            self.x0[ix_['x177']] = float(1251.6150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x177']),float(1251.6150000)))
        if('y177' in ix_):
            self.x0[ix_['y177']] = float(62993.389650)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y177']),float(62993.389650)))
        if('w178' in ix_):
            self.x0[ix_['w178']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w178']),float(13.000000000)))
        if('sigt178' in ix_):
            self.x0[ix_['sigt178']] = float(15.036517615)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt178']),float(15.036517615)))
        if('sigr178' in ix_):
            self.x0[ix_['sigr178']] = float(19.101702555)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr178']),float(19.101702555)))
        if('x178' in ix_):
            self.x0[ix_['x178']] = float(1258.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x178']),float(1258.1150000)))
        if('y178' in ix_):
            self.x0[ix_['y178']] = float(63093.229145)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y178']),float(63093.229145)))
        if('w179' in ix_):
            self.x0[ix_['w179']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w179']),float(13.000000000)))
        if('sigt179' in ix_):
            self.x0[ix_['sigt179']] = float(14.383974972)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt179']),float(14.383974972)))
        if('sigr179' in ix_):
            self.x0[ix_['sigr179']] = float(18.801977014)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr179']),float(18.801977014)))
        if('x179' in ix_):
            self.x0[ix_['x179']] = float(1264.6150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x179']),float(1264.6150000)))
        if('y179' in ix_):
            self.x0[ix_['y179']] = float(63188.845746)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y179']),float(63188.845746)))
        if('w180' in ix_):
            self.x0[ix_['w180']] = float(13.000000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['w180']),float(13.000000000)))
        if('sigr180' in ix_):
            self.x0[ix_['sigr180']] = float(18.500000000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigr180']),float(18.500000000)))
        if('sigt180' in ix_):
            self.x0[ix_['sigt180']] = float(13.725650411)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['sigt180']),float(13.725650411)))
        if('x180' in ix_):
            self.x0[ix_['x180']] = float(1271.1150000)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['x180']),float(1271.1150000)))
        if('y180' in ix_):
            self.x0[ix_['y180']] = float(63280.202028)
        else:
            self.y0  = (
                  arrset(self.y0,findfirst(self.congrps,lambda x:x==ig_['y180']),float(63280.202028)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'en2PR', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for k in range(int(v_['0']),int(v_['K'])+1):
            ename = 'WSR'+str(k)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_["en2PR"])
            vname = 'w'+str(k)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'sigr'+str(k)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for k in range(int(v_['0']),int(v_['K'])+1):
            ename = 'WST'+str(k)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en2PR')
            ielftype = arrset(ielftype,ie,iet_["en2PR"])
            vname = 'w'+str(k)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'sigt'+str(k)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        v_['rk'] = v_['ri']
        for k in range(int(v_['0']),int(v_['K-1'])+1):
            v_['rk+1'] = v_['rk']+v_['dr']
            v_['k+1'] = 1+k
            v_['-rk'] = -1.0*v_['rk']
            ig = ig_['SR'+str(k)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WSR'+str(k)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-rk']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WSR'+str(int(v_['k+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['rk+1']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WST'+str(k)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-dr/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WST'+str(int(v_['k+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-dr/2']))
            ig = ig_['STAy'+str(k)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WST'+str(k)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-dr/2']))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['WST'+str(int(v_['k+1']))])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(v_['-dr/2']))
            v_['rk'] = v_['rk+1']
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        self.objlower = 5.0
# LO SOLUTION            7.872067544
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
        self.pbclass   = "C-CLQR2-RN-905-1081"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def en2PR(self, nargout,*args):

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

