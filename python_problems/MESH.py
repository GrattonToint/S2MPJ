from s2mpjlib import *
class  MESH(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The goodness of a finite element grid is characterized by the
#    smallest angle of all triangles. Given a triangulation of a domain 
#    in R**2. Find a topological equivalent triangulation so, that the 
#    smallest angel becomes as large as possible. Topological equivalent 
#    means shifting the edges of the grid only in such a way that 
#    neighbouring triangles remain neighbours. 
# 
#    Source: Prof. Dr. Michael Kraetzschmar, Institut fuer Angewandte 
#            Mathematik der Fachhochschule Flensburg, Kanzleistrasse 91-93, 
#            D-24943 FLENSBURG, GERMANY
# 
#    SIF input: Prof. Dr. Michael Kraetzschmar
# 
#    classification = "C-COOR2-AY-41-48"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MESH'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['omega1'] = 1.0e3
        v_['omega2'] = -1.0e3
        v_['omega3'] = -1.0e5
        v_['s'] = 0.700000
        v_['pi'] = np.arccos(-1.0)
        v_['sqrt3/2'] = np.sqrt(0.75)
        v_['h'] = v_['sqrt3/2']*v_['s']
        v_['drei'] = 3.0
        v_['pi/3'] = v_['pi']/v_['drei']
        v_['1'] = 1
        v_['np'] = 5
        v_['nk'] = 8
        v_['nd'] = 4
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for i in range(int(v_['1']),int(v_['np'])+1):
            [iv,ix_,_] = s2mpj_ii('x'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'x'+str(i))
            [iv,ix_,_] = s2mpj_ii('y'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'y'+str(i))
        for i in range(int(v_['1']),int(v_['nk'])+1):
            [iv,ix_,_] = s2mpj_ii('l'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'l'+str(i))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [iv,ix_,_] = s2mpj_ii('alpha'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'alpha'+str(i))
            [iv,ix_,_] = s2mpj_ii('beta'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'beta'+str(i))
            [iv,ix_,_] = s2mpj_ii('gamma'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'gamma'+str(i))
            [iv,ix_,_] = s2mpj_ii('delta'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'delta'+str(i))
            [iv,ix_,_] = s2mpj_ii('f'+str(i),ix_)
            self.xnames=arrset(self.xnames,iv,'f'+str(i))
        [iv,ix_,_] = s2mpj_ii('deltamin',ix_)
        self.xnames=arrset(self.xnames,iv,'deltamin')
        [iv,ix_,_] = s2mpj_ii('fmin',ix_)
        self.xnames=arrset(self.xnames,iv,'fmin')
        [iv,ix_,_] = s2mpj_ii('fmax',ix_)
        self.xnames=arrset(self.xnames,iv,'fmax')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('obj1',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['deltamin']])
        valA = np.append(valA,float(1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['omega1']))
        [ig,ig_,_] = s2mpj_ii('obj2',ig_)
        gtype = arrset(gtype,ig,'<>')
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['fmax']])
        valA = np.append(valA,float(1.0))
        irA  = np.append(irA,[ig])
        icA  = np.append(icA,[ix_['fmin']])
        valA = np.append(valA,float(-1.0))
        self.gscale = arrset(self.gscale,ig,float(v_['omega2']))
        [ig,ig_,_] = s2mpj_ii('obj3',ig_)
        gtype = arrset(gtype,ig,'<>')
        self.gscale = arrset(self.gscale,ig,float(v_['omega3']))
        for i in range(int(v_['1']),int(v_['nk'])+1):
            [ig,ig_,_] = s2mpj_ii('seit'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'seit'+str(i))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [ig,ig_,_] = s2mpj_ii('skal'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'skal'+str(i))
            [ig,ig_,_] = s2mpj_ii('skbe'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'skbe'+str(i))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [ig,ig_,_] = s2mpj_ii('doppf'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'doppf'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['f'+str(i)]])
            valA = np.append(valA,float(-1.0))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [ig,ig_,_] = s2mpj_ii('wisum'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'wisum'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['alpha'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['beta'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['gamma'+str(i)]])
            valA = np.append(valA,float(1.0))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [ig,ig_,_] = s2mpj_ii('alphd'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'alphd'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['alpha'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['delta'+str(i)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('betad'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'betad'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['beta'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['delta'+str(i)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('gammd'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'gammd'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['gamma'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['delta'+str(i)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('deltd'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'deltd'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['delta'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['deltamin']])
            valA = np.append(valA,float(-1.0))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            [ig,ig_,_] = s2mpj_ii('fmind'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'fmind'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['f'+str(i)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['fmin']])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('fmaxd'+str(i),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'fmaxd'+str(i))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['fmax']])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['f'+str(i)]])
            valA = np.append(valA,float(-1.0))
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
        for i in range(int(v_['1']),int(v_['nd'])+1):
            self.gconst = arrset(self.gconst,ig_['wisum'+str(i)],float(v_['pi']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xupper[ix_['deltamin']] = v_['pi']
        self.xlower[ix_['x1']] = 0.000000
        self.xupper[ix_['x1']] = 0.000000
        self.xlower[ix_['y1']] = 0.000000
        self.xupper[ix_['y1']] = 0.000000
        self.xlower[ix_['x2']] = 0.000000
        self.xupper[ix_['x2']] = 0.000000
        self.xlower[ix_['y2']] = 1.000000
        self.xupper[ix_['y2']] = 1.000000
        self.xlower[ix_['x3']] = 1.000000
        self.xupper[ix_['x3']] = 1.000000
        self.xlower[ix_['y3']] = 1.000000
        self.xupper[ix_['y3']] = 1.000000
        self.xlower[ix_['x4']] = 1.000000
        self.xupper[ix_['x4']] = 1.000000
        self.xlower[ix_['y4']] = 0.000000
        self.xupper[ix_['y4']] = 0.000000
        self.xlower[ix_['x5']] = -float('Inf')
        self.xupper[ix_['x5']] = +float('Inf')
        self.xlower[ix_['y5']] = -float('Inf')
        self.xupper[ix_['y5']] = +float('Inf')
        for i in range(int(v_['1']),int(v_['nd'])+1):
            self.xupper[ix_['alpha'+str(i)]] = v_['pi']
            self.xupper[ix_['beta'+str(i)]] = v_['pi']
            self.xupper[ix_['gamma'+str(i)]] = v_['pi']
            self.xupper[ix_['delta'+str(i)]] = v_['pi']
        self.xupper[ix_['deltamin']] = v_['pi']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['x1']] = float(0.000000)
        self.x0[ix_['y1']] = float(0.000000)
        self.x0[ix_['x2']] = float(0.000000)
        self.x0[ix_['y2']] = float(1.000000)
        self.x0[ix_['x3']] = float(1.000000)
        self.x0[ix_['y3']] = float(1.000000)
        self.x0[ix_['x4']] = float(1.000000)
        self.x0[ix_['y4']] = float(0.000000)
        self.x0[ix_['x5']] = float(0.350000)
        self.x0[ix_['y5']] = float(0.606218)
        self.x0[ix_['l1']] = float(0.700000)
        self.x0[ix_['l2']] = float(0.526844)
        self.x0[ix_['l3']] = float(1.000000)
        self.x0[ix_['l4']] = float(0.759977)
        self.x0[ix_['l5']] = float(1.000000)
        self.x0[ix_['l6']] = float(0.888819)
        self.x0[ix_['l7']] = float(1.000000)
        self.x0[ix_['l8']] = float(1.000000)
        self.x0[ix_['alpha1']] = float(0.523599)
        self.x0[ix_['beta1']] = float(1.891392)
        self.x0[ix_['gamma1']] = float(0.726602)
        self.x0[ix_['delta1']] = float(0.523599)
        self.x0[ix_['f1']] = float(0.350000)
        self.x0[ix_['alpha2']] = float(0.844195)
        self.x0[ix_['beta2']] = float(1.752711)
        self.x0[ix_['gamma2']] = float(0.544687)
        self.x0[ix_['delta2']] = float(0.544687)
        self.x0[ix_['f2']] = float(0.393782)
        self.x0[ix_['alpha3']] = float(1.026109)
        self.x0[ix_['beta3']] = float(1.295247)
        self.x0[ix_['gamma3']] = float(0.820236)
        self.x0[ix_['delta3']] = float(0.820236)
        self.x0[ix_['f3']] = float(0.650000)
        self.x0[ix_['alpha4']] = float(0.750560)
        self.x0[ix_['beta4']] = float(1.343835)
        self.x0[ix_['gamma4']] = float(1.047198)
        self.x0[ix_['delta4']] = float(0.750560)
        self.x0[ix_['f4']] = float(0.606218)
        self.x0[ix_['deltamin']] = float(0.523599)
        self.x0[ix_['fmin']] = float(0.350000)
        self.x0[ix_['fmax']] = float(0.650000)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ediffsq', iet_)
        elftv = loaset(elftv,it,0,'winkel')
        elftv = loaset(elftv,it,1,'minwi')
        [it,iet_,_] = s2mpj_ii( 'elaenge', iet_)
        elftv = loaset(elftv,it,0,'lp1x')
        elftv = loaset(elftv,it,1,'lp1y')
        elftv = loaset(elftv,it,2,'lp2x')
        elftv = loaset(elftv,it,3,'lp2y')
        elftv = loaset(elftv,it,4,'llaenge')
        [it,iet_,_] = s2mpj_ii( 'evekprod', iet_)
        elftv = loaset(elftv,it,0,'vp1x')
        elftv = loaset(elftv,it,1,'vp1y')
        elftv = loaset(elftv,it,2,'vp2x')
        elftv = loaset(elftv,it,3,'vp2y')
        elftv = loaset(elftv,it,4,'vp3x')
        elftv = loaset(elftv,it,5,'vp3y')
        [it,iet_,_] = s2mpj_ii( 'esklprod', iet_)
        elftv = loaset(elftv,it,0,'sp1x')
        elftv = loaset(elftv,it,1,'sp1y')
        elftv = loaset(elftv,it,2,'sp2x')
        elftv = loaset(elftv,it,3,'sp2y')
        elftv = loaset(elftv,it,4,'sp3x')
        elftv = loaset(elftv,it,5,'sp3y')
        [it,iet_,_] = s2mpj_ii( 'ecosprod', iet_)
        elftv = loaset(elftv,it,0,'cl1')
        elftv = loaset(elftv,it,1,'cl2')
        elftv = loaset(elftv,it,2,'cw')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'aldsq'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ediffsq')
            ielftype = arrset(ielftype,ie,iet_["ediffsq"])
            vname = 'alpha'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='winkel')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'delta'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='minwi')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'bedsq'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ediffsq')
            ielftype = arrset(ielftype,ie,iet_["ediffsq"])
            vname = 'beta'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='winkel')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'delta'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='minwi')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'gadsq'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ediffsq')
            ielftype = arrset(ielftype,ie,iet_["ediffsq"])
            vname = 'gamma'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='winkel')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'delta'+str(i)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='minwi')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nk'])+1):
            ename = 'laeng'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'elaenge')
            ielftype = arrset(ielftype,ie,iet_["elaenge"])
        ename = 'laeng1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng5'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng6'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng7'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'laeng8'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='lp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='llaenge')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'sal'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'esklprod')
            ielftype = arrset(ielftype,ie,iet_["esklprod"])
        ename = 'sal1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sal2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sal3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sal4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'cal'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ecosprod')
            ielftype = arrset(ielftype,ie,iet_["ecosprod"])
        ename = 'cal1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'alpha1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cal2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'alpha2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cal3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'alpha3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cal4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'alpha4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'sbe'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'esklprod')
            ielftype = arrset(ielftype,ie,iet_["esklprod"])
        ename = 'sbe1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sbe2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sbe3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'sbe4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='sp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'cbe'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ecosprod')
            ielftype = arrset(ielftype,ie,iet_["ecosprod"])
        ename = 'cbe1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'beta1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cbe2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'beta2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cbe3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'beta3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'cbe4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'l6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'l1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cl2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'beta4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='cw')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ename = 'flae'+str(i)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'evekprod')
            ielftype = arrset(ielftype,ie,iet_["evekprod"])
        ename = 'flae1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'flae2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'flae3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'flae4'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        vname = 'x4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp1y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp2y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'x1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3x')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'y1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='vp3y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gsquare',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['obj1']
        self.grftype = arrset(self.grftype,ig,'gsquare')
        ig = ig_['obj2']
        self.grftype = arrset(self.grftype,ig,'gsquare')
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ig = ig_['obj3']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['aldsq'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['bedsq'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['gadsq'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        for i in range(int(v_['1']),int(v_['nk'])+1):
            ig = ig_['seit'+str(i)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['laeng'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ig = ig_['skal'+str(i)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['sal'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['cal'+str(i)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
            ig = ig_['skbe'+str(i)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['sbe'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
            posel = posel+1
            self.grelt = loaset(self.grelt,ig,posel,ie_['cbe'+str(i)])
            self.grelw = loaset(self.grelw,ig,posel,float(-1.0))
        for i in range(int(v_['1']),int(v_['nd'])+1):
            ig = ig_['doppf'+str(i)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['flae'+str(i)])
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              5.9213448D-4
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.clower[np.arange(self.nle+self.neq,self.m)] = np.zeros((self.nge,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-COOR2-AY-41-48"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ediffsq(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        f_   = IV_[0]*IV_[0]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[0]+IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def elaenge(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((3,5))
        IV_ = np.zeros(3)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,4] = U_[2,4]+1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        f_   = IV_[2]*IV_[2]-IV_[0]*IV_[0]-IV_[1]*IV_[1]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[2] = IV_[2]+IV_[2]
            g_[0] = -(IV_[0]+IV_[0])
            g_[1] = -(IV_[1]+IV_[1])
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[2,2] = 2.0
                H_[2,0] = 0.0
                H_[0,2] = H_[2,0]
                H_[2,1] = 0.0
                H_[1,2] = H_[2,1]
                H_[0,1] = 0.0
                H_[1,0] = H_[0,1]
                H_[0,0] = -2.0
                H_[1,1] = -2.0
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def evekprod(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,6))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,4] = U_[2,4]+1
        U_[2,2] = U_[2,2]-1
        U_[3,5] = U_[3,5]+1
        U_[3,3] = U_[3,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        f_   = IV_[2]*IV_[1]-IV_[0]*IV_[3]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -IV_[3]
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_[3] = -IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[2,1] = 1.0
                H_[1,2] = H_[2,1]
                H_[0,3] = -1.0
                H_[3,0] = H_[0,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def esklprod(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((4,6))
        IV_ = np.zeros(4)
        U_[0,0] = U_[0,0]+1
        U_[0,2] = U_[0,2]-1
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,4] = U_[2,4]+1
        U_[2,2] = U_[2,2]-1
        U_[3,5] = U_[3,5]+1
        U_[3,3] = U_[3,3]-1
        IV_[0] = to_scalar(U_[0:1,:].dot(EV_))
        IV_[1] = to_scalar(U_[1:2,:].dot(EV_))
        IV_[2] = to_scalar(U_[2:3,:].dot(EV_))
        IV_[3] = to_scalar(U_[3:4,:].dot(EV_))
        f_   = IV_[0]*IV_[2]+IV_[1]*IV_[3]
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[2]
            g_[1] = IV_[3]
            g_[2] = IV_[0]
            g_[3] = IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,2] = 1.0
                H_[2,0] = H_[0,2]
                H_[1,3] = 1.0
                H_[3,1] = H_[1,3]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ecosprod(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        cosa = np.cos(EV_[2,0])
        sina = np.sin(EV_[2,0])
        prod2 = EV_[0,0]*EV_[1,0]
        prod3 = prod2*cosa
        f_   = prod3
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1,0]*cosa
            g_[1] = EV_[0,0]*cosa
            g_[2] = -sina*prod2
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = cosa
                H_[1,0] = H_[0,1]
                H_[0,2] = -sina*EV_[1,0]
                H_[2,0] = H_[0,2]
                H_[1,2] = -sina*EV_[0,0]
                H_[2,1] = H_[1,2]
                H_[2,2] = -prod3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gsquare(self,nargout,*args):

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

