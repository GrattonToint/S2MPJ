from s2xlib import *
class  BAmL1(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BAmL1
#    *********
# 
#    Bundle Adjustment problem from reconstructive geometry in which
#    a collection of photographs is used to determine the position of
#    a set of observed points. Each observed point is seen via its
#    two-dimensional projections on a subset of the photographs. The
#    solution is found by solvng a large nonlinear least-squares problem.
#    This variant is given as an inconsistent set of nonlinear equations.
# 
#    Source: data from the Bundle Adjustment in the Large
#    project, http://grail.cs.washington.edu/projects/bal/
# 
#    Ladybug datasets (single image extracted)
# 
#    SIF input: Nick Gould, November 2016
# 
#    classification = "NOR2-MN-57-12"
# 
#    Number of images
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BAmL1'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'BAmL1'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['nuimages'] = 49
        v_['nupoints'] = 1
        v_['nuobservs'] = 6
        v_['1'] = 1
        v_['O1'] = 1.0
        v_['O2'] = 2.0
        v_['O3'] = 4.0
        v_['O4'] = 27.0
        v_['O5'] = 30.0
        v_['O6'] = 37.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['nupoints'])+1):
            [iv,ix_,_] = s2x_ii('X'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'X'+str(I))
            [iv,ix_,_] = s2x_ii('Y'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Y'+str(I))
            [iv,ix_,_] = s2x_ii('Z'+str(I),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Z'+str(I))
        for J in range(int(v_['1']),int(v_['nuobservs'])+1):
            v_['RI'] = v_['O'+str(J)]
            v_['I'] = int(np.fix(v_['RI']))
            [iv,ix_,_] = s2x_ii('RX'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'RX'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('RY'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'RY'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('RZ'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'RZ'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('TX'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'TX'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('TY'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'TY'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('TZ'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'TZ'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('KA'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'KA'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('KB'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'KB'+str(int(v_['I'])))
            [iv,ix_,_] = s2x_ii('F'+str(int(v_['I'])),ix_)
            pb.xnames=arrset(pb.xnames,iv,'F'+str(int(v_['I'])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['nuobservs'])+1):
            [ig,ig_,_] = s2x_ii('RX'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'RX'+str(I))
            [ig,ig_,_] = s2x_ii('RY'+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'RY'+str(I))
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
        pbm.gconst = arrset(pbm.gconst,ig_['RX1'],float(-332.65))
        pbm.gconst = arrset(pbm.gconst,ig_['RY1'],float(262.09))
        pbm.gconst = arrset(pbm.gconst,ig_['RX2'],float(-199.76))
        pbm.gconst = arrset(pbm.gconst,ig_['RY2'],float(166.7))
        pbm.gconst = arrset(pbm.gconst,ig_['RX3'],float(-253.06))
        pbm.gconst = arrset(pbm.gconst,ig_['RY3'],float(202.27))
        pbm.gconst = arrset(pbm.gconst,ig_['RX4'],float(58.13))
        pbm.gconst = arrset(pbm.gconst,ig_['RY4'],float(271.89))
        pbm.gconst = arrset(pbm.gconst,ig_['RX5'],float(238.22))
        pbm.gconst = arrset(pbm.gconst,ig_['RY5'],float(237.37))
        pbm.gconst = arrset(pbm.gconst,ig_['RX6'],float(317.55))
        pbm.gconst = arrset(pbm.gconst,ig_['RY6'],float(221.15))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('X1' in ix_):
            pb.x0[ix_['X1']] = float(-.6120001572)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X1']),float(-.6120001572)))
        if('Y1' in ix_):
            pb.x0[ix_['Y1']] = float(.57175904776)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Y1']),float(.57175904776)))
        if('Z1' in ix_):
            pb.x0[ix_['Z1']] = float(-1.847081276)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Z1']),float(-1.847081276)))
        if('RX1' in ix_):
            pb.x0[ix_['RX1']] = float(.01574151594)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX1']),float(.01574151594)))
        if('RY1' in ix_):
            pb.x0[ix_['RY1']] = float(-.0127909362)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY1']),float(-.0127909362)))
        if('RZ1' in ix_):
            pb.x0[ix_['RZ1']] = float(-.0044008498)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ1']),float(-.0044008498)))
        if('TX1' in ix_):
            pb.x0[ix_['TX1']] = float(-.0340938396)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX1']),float(-.0340938396)))
        if('TY1' in ix_):
            pb.x0[ix_['TY1']] = float(-.107513871)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY1']),float(-.107513871)))
        if('TZ1' in ix_):
            pb.x0[ix_['TZ1']] = float(1.1202240291)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ1']),float(1.1202240291)))
        if('KA1' in ix_):
            pb.x0[ix_['KA1']] = float(-3.177064E-7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA1']),float(-3.177064E-7)))
        if('KB1' in ix_):
            pb.x0[ix_['KB1']] = float(5.882049E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB1']),float(5.882049E-13)))
        if('F1' in ix_):
            pb.x0[ix_['F1']] = float(399.75152639)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F1']),float(399.75152639)))
        if('RX2' in ix_):
            pb.x0[ix_['RX2']] = float(.01597732412)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX2']),float(.01597732412)))
        if('RY2' in ix_):
            pb.x0[ix_['RY2']] = float(-.0252244646)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY2']),float(-.0252244646)))
        if('RZ2' in ix_):
            pb.x0[ix_['RZ2']] = float(-.0094001416)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ2']),float(-.0094001416)))
        if('TX2' in ix_):
            pb.x0[ix_['TX2']] = float(-.0085667661)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX2']),float(-.0085667661)))
        if('TY2' in ix_):
            pb.x0[ix_['TY2']] = float(-.1218804907)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY2']),float(-.1218804907)))
        if('TZ2' in ix_):
            pb.x0[ix_['TZ2']] = float(.7190133075)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ2']),float(.7190133075)))
        if('KA2' in ix_):
            pb.x0[ix_['KA2']] = float(-3.780477E-7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA2']),float(-3.780477E-7)))
        if('KB2' in ix_):
            pb.x0[ix_['KB2']] = float(9.307431E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB2']),float(9.307431E-13)))
        if('F2' in ix_):
            pb.x0[ix_['F2']] = float(402.01753386)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F2']),float(402.01753386)))
        if('RX4' in ix_):
            pb.x0[ix_['RX4']] = float(.01484625118)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX4']),float(.01484625118)))
        if('RY4' in ix_):
            pb.x0[ix_['RY4']] = float(-.0210628994)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY4']),float(-.0210628994)))
        if('RZ4' in ix_):
            pb.x0[ix_['RZ4']] = float(-.001166948)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ4']),float(-.001166948)))
        if('TX4' in ix_):
            pb.x0[ix_['TX4']] = float(-.0249509707)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX4']),float(-.0249509707)))
        if('TY4' in ix_):
            pb.x0[ix_['TY4']] = float(-.1139847055)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY4']),float(-.1139847055)))
        if('TZ4' in ix_):
            pb.x0[ix_['TZ4']] = float(.92166020737)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ4']),float(.92166020737)))
        if('KA4' in ix_):
            pb.x0[ix_['KA4']] = float(-3.295265E-7)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA4']),float(-3.295265E-7)))
        if('KB4' in ix_):
            pb.x0[ix_['KB4']] = float(6.732885E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB4']),float(6.732885E-13)))
        if('F4' in ix_):
            pb.x0[ix_['F4']] = float(400.40175368)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F4']),float(400.40175368)))
        if('RX27' in ix_):
            pb.x0[ix_['RX27']] = float(.01991666998)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX27']),float(.01991666998)))
        if('RY27' in ix_):
            pb.x0[ix_['RY27']] = float(-1.22433082)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY27']),float(-1.22433082)))
        if('RZ27' in ix_):
            pb.x0[ix_['RZ27']] = float(.0119988756)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ27']),float(.0119988756)))
        if('TX27' in ix_):
            pb.x0[ix_['TX27']] = float(-1.411897512)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX27']),float(-1.411897512)))
        if('TY27' in ix_):
            pb.x0[ix_['TY27']] = float(-.1148065151)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY27']),float(-.1148065151)))
        if('TZ27' in ix_):
            pb.x0[ix_['TZ27']] = float(.44915582738)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ27']),float(.44915582738)))
        if('KA27' in ix_):
            pb.x0[ix_['KA27']] = float(5.95875E-8)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA27']),float(5.95875E-8)))
        if('KB27' in ix_):
            pb.x0[ix_['KB27']] = float(-2.48391E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB27']),float(-2.48391E-13)))
        if('F27' in ix_):
            pb.x0[ix_['F27']] = float(407.03024568)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F27']),float(407.03024568)))
        if('RX30' in ix_):
            pb.x0[ix_['RX30']] = float(.02082242153)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX30']),float(.02082242153)))
        if('RY30' in ix_):
            pb.x0[ix_['RY30']] = float(-1.238434791)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY30']),float(-1.238434791)))
        if('RZ30' in ix_):
            pb.x0[ix_['RZ30']] = float(.01389314763)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ30']),float(.01389314763)))
        if('TX30' in ix_):
            pb.x0[ix_['TX30']] = float(-1.049686225)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX30']),float(-1.049686225)))
        if('TY30' in ix_):
            pb.x0[ix_['TY30']] = float(-.1299513286)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY30']),float(-.1299513286)))
        if('TZ30' in ix_):
            pb.x0[ix_['TZ30']] = float(.33798380231)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ30']),float(.33798380231)))
        if('KA30' in ix_):
            pb.x0[ix_['KA30']] = float(4.5673127E-8)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA30']),float(4.5673127E-8)))
        if('KB30' in ix_):
            pb.x0[ix_['KB30']] = float(-1.79243E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB30']),float(-1.79243E-13)))
        if('F30' in ix_):
            pb.x0[ix_['F30']] = float(405.91764962)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F30']),float(405.91764962)))
        if('RX37' in ix_):
            pb.x0[ix_['RX37']] = float(.01658816461)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RX37']),float(.01658816461)))
        if('RY37' in ix_):
            pb.x0[ix_['RY37']] = float(-1.247226838)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RY37']),float(-1.247226838)))
        if('RZ37' in ix_):
            pb.x0[ix_['RZ37']] = float(.01846788123)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['RZ37']),float(.01846788123)))
        if('TX37' in ix_):
            pb.x0[ix_['TX37']] = float(-.8617315756)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TX37']),float(-.8617315756)))
        if('TY37' in ix_):
            pb.x0[ix_['TY37']] = float(-.1321089362)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TY37']),float(-.1321089362)))
        if('TZ37' in ix_):
            pb.x0[ix_['TZ37']] = float(.28256800868)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['TZ37']),float(.28256800868)))
        if('KA37' in ix_):
            pb.x0[ix_['KA37']] = float(4.7465711E-8)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KA37']),float(4.7465711E-8)))
        if('KB37' in ix_):
            pb.x0[ix_['KB37']] = float(-1.50881E-13)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['KB37']),float(-1.50881E-13)))
        if('F37' in ix_):
            pb.x0[ix_['F37']] = float(404.73590637)
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['F37']),float(404.73590637)))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eE', iet_)
        elftv = loaset(elftv,it,0,'RX')
        elftv = loaset(elftv,it,1,'RY')
        elftv = loaset(elftv,it,2,'RZ')
        elftv = loaset(elftv,it,3,'X')
        elftv = loaset(elftv,it,4,'Y')
        elftv = loaset(elftv,it,5,'Z')
        elftv = loaset(elftv,it,6,'TX')
        elftv = loaset(elftv,it,7,'TY')
        elftv = loaset(elftv,it,8,'TZ')
        elftv = loaset(elftv,it,9,'KA')
        elftv = loaset(elftv,it,10,'KB')
        elftv = loaset(elftv,it,11,'F')
        elftp = []
        elftp = loaset(elftp,it,0,'YRES')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        ename = 'EX1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY1'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EX2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY2'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F2'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EX3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY3'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F4'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EX4'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY4'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F27'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EX5'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY5'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F30'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        ename = 'EX6'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(0.0))
        ename = 'EY6'
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'eE')
        ielftype = arrset(ielftype, ie, iet_["eE"])
        vname = 'X1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Y1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'Z1'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Z')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RX37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RY37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'RZ37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='RZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TX37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TX')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TY37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TY')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'TZ37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='TZ')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KA37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KA')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'KB37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='KB')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        vname = 'F37'
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='F')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        posep = find(elftp[ielftype[ie]],lambda x:x=='YRES')
        pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['nuobservs'])+1):
            ig = ig_['RX'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EX'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            ig = ig_['RY'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['EY'+str(I)])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NOR2-MN-57-12"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
