from s2mpjlib import *
class  SANTA(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SANTA
#    --------
# 
#    The Santa problem as suggested in a Christmas competition
#    by Jens Jensen (Scientific Computing, STFC). To quote Jens,
# 
#    Santa and His Elves
# 
#    SCD Christmas programming challenge 2016
# 
#    Christmas has come to the Santa Claus Department – or rather, the SCD
#    is coming to Christmas. Santa is flying around the world, presently
#    presenting presents. Ho, ho, ho! No striking air crew on Santa’s sleigh!
#    No airport strikes on the North Pole.
# 
#    For the purpose of this exercise, the Earth is round as a perfect ball,
#    with radius precisely 6,371,000 metres. However, everything is at the
#    same longitude and latitude as the “real” Earth. So for example, the
#    Greenwich observatory is at 51°28'40"N 0°00'04"W both on the “real”
#    Earth and on Santa’s Earth. (Also ignore rotation of the Earth and
#    anything practical like that.)
# 
#    Santa sets off from the North Pole along 2°6'57.6" E bearing south
#    (obviously), and also bearing presents (obviously). Whenever Santa
#    leaves a location, he leaves an elf behind, in order to help unwrapping
#    presents; the elf does this and then flies out independently to meet up
#    with Santa at the next location - this means that Santa only needs two
#    Elves. Here’s how:
# 
#    1. Santa leaves the North Pole, setting out for location A. Elf 1 is
#    left behind (in this particular case, not to unwrap presents, but to
#    turn the lights off, and ensure the oven is off – it's elf'n'safety,
#    you know.)
# 
#    2. Santa arrives in location A and hands them their present. Now Elf 2
#    is with Santa; Elf 1 is still at the NP.
# 
#    3. Santa leaves location A, leaving behind Elf 2. Santa flies on to
#    location B; Elf 1, who remained at the North Pole, also flies to B and
#    meets Santa there; Elf 2 is left behind at A.
# 
#    4. Santa arrives at location B along with Elf 1, and hands out
#    presents. Santa then leaves location B to fly to C, leaving behind Elf 1
#    at location B. Meanwhile Elf 2, having finished helping at location A,
#    leaves location A to fly on to C, to meet Santa there.
# 
#    5. Santa arrives from B at location C; Elf 2 also arrives into C from
#    location A. Elf 1 remains at B until Santa flies onward to location D.
# 
#    6. At the last hop, Santa needs a rest and flies to 31°46'42.4" S
#    144°46'12.9" W.  The Elves also fly to this location - maps show no land
#    here but it is hidden. Either that or we got the coordinates wrong.
#    In either case Santa and elves do fly to this location.
# 
#    The following table shows the distance of Santa's hops, as well as those
#    of the elves, with the distance given in metres:
# 
#    Who     Hop  Distance travelled
#    Santa   1    5405238
#            2    623852
#            3    1005461
#            4    7470967
#            5    3632559
#            6    10206818
#            7    7967212
#            8    5896361
#            9    8337266
#            10   13019505
#            11   8690818
#            12   8971302
#    Elf1    1    4866724
#            2    6833740
#            3    13489586
#            4    9195575
#            5    9704793
#            6    12498127
#    Elf2    1    1375828
#            2    4917407
#            3    10617953
#            4    10996150
#            5    7901038
#            6    8971302
# 
#    What is Santa’s route?  What sort of presents is he carrying?
# 
#    Bonus question: did you really need to know the starting direction?
# 
#    Added by Nick: the problem has many local minimizers of the sum of squares
#    of infeasibility, but it is only the solution with zero residuals
#    that is of interest.
# 
#    SIF input: Nick Gould, Dec 2016.
# 
#    classification = "C-CNOR2-AN-21-23"
# 
#    Number of stops on Santa's path (path goes from index 0 to 12)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SANTA'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['S'] = 12
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['180.0'] = 180.0
        v_['S-1'] = -1+v_['S']
        v_['RS'] = float(v_['S'])
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['PI/180'] = v_['PI']/v_['180.0']
        v_['PHI0'] = 90.0
        v_['LAM0'] = 0.0
        v_['PHI12'] = -31.77844444
        v_['LAM12'] = -144.77025
        v_['LAM1'] = 2.116
        v_['PHI0'] = v_['PHI0']*v_['PI/180']
        v_['LAM0'] = v_['LAM0']*v_['PI/180']
        v_['PHI12'] = v_['PHI12']*v_['PI/180']
        v_['LAM12'] = v_['LAM12']*v_['PI/180']
        v_['LAM1'] = v_['LAM1']*v_['PI/180']
        v_['DPHI'] = v_['PHI12']-v_['PHI0']
        v_['DLAM'] = v_['LAM12']-v_['LAM0']
        v_['DPHI/S'] = v_['DPHI']/v_['RS']
        v_['DLAM/S'] = v_['DLAM']/v_['RS']
        v_['RADIUS'] = 6371000.0
        v_['D0,1'] = 5405238.0
        v_['D0,2'] = 4866724.0
        v_['D1,2'] = 623852.0
        v_['D1,3'] = 1375828.0
        v_['D2,3'] = 1005461.0
        v_['D2,4'] = 6833740.0
        v_['D3,4'] = 7470967.0
        v_['D3,5'] = 4917407.0
        v_['D4,5'] = 3632559.0
        v_['D4,6'] = 13489586.0
        v_['D5,6'] = 10206818.0
        v_['D5,7'] = 10617953.0
        v_['D6,7'] = 7967212.0
        v_['D6,8'] = 9195575.0
        v_['D7,8'] = 5896361.0
        v_['D7,9'] = 10996150.0
        v_['D8,9'] = 8337266.0
        v_['D8,10'] = 9704793.0
        v_['D9,10'] = 13019505.0
        v_['D9,11'] = 7901038.0
        v_['D10,11'] = 8690818.0
        v_['D10,12'] = 12498127.0
        v_['D11,12'] = 8971302.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('PHI1',ix_)
        self.xnames=arrset(self.xnames,iv,'PHI1')
        for I in range(int(v_['2']),int(v_['S-1'])+1):
            [iv,ix_,_] = s2mpj_ii('PHI'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'PHI'+str(I))
            [iv,ix_,_] = s2mpj_ii('LAM'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'LAM'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('R0,1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'R0,1')
        for I in range(int(v_['2']),int(v_['S'])+1):
            v_['I1'] = -1+I
            v_['I2'] = -2+I
            [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['I2']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(int(v_['I2']))+','+str(I))
            [ig,ig_,_] = s2mpj_ii('R'+str(int(v_['I1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'R'+str(int(v_['I1']))+','+str(I))
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
        v_['D/RAD'] = v_['D0,1']/v_['RADIUS']
        v_['CD/RAD'] = np.cos(v_['D/RAD'])
        self.gconst = arrset(self.gconst,ig_['R0,1'],float(v_['CD/RAD']))
        for I in range(int(v_['2']),int(v_['S'])+1):
            v_['I2'] = -2+I
            v_['D/RAD'] = v_['D'+str(int(v_['I2']))+','+str(I)]/v_['RADIUS']
            v_['CD/RAD'] = np.cos(v_['D/RAD'])
            self.gconst  = (
                  arrset(self.gconst,ig_['R'+str(int(v_['I2']))+','+str(I)],float(v_['CD/RAD'])))
            v_['I1'] = -1+I
            v_['D/RAD'] = v_['D'+str(int(v_['I1']))+','+str(I)]/v_['RADIUS']
            v_['CD/RAD'] = np.cos(v_['D/RAD'])
            self.gconst  = (
                  arrset(self.gconst,ig_['R'+str(int(v_['I1']))+','+str(I)],float(v_['CD/RAD'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-1000.0)
        self.xupper = np.full((self.n,1),1000.0)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['PHI1']] = float(0.7223835215)
        self.x0[ix_['PHI2']] = float(0.8069093428)
        self.x0[ix_['LAM2']] = float(-0.031657133)
        self.x0[ix_['PHI3']] = float(0.9310164154)
        self.x0[ix_['LAM3']] = float(0.1199353230)
        self.x0[ix_['PHI4']] = float(6.6067392710)
        self.x0[ix_['LAM4']] = float(-1.214314477)
        self.x0[ix_['PHI5']] = float(-3.530946794)
        self.x0[ix_['LAM5']] = float(2.5329493980)
        self.x0[ix_['PHI6']] = float(-9.798251905)
        self.x0[ix_['LAM6']] = float(4.3021328700)
        self.x0[ix_['PHI7']] = float(14.632267534)
        self.x0[ix_['LAM7']] = float(-12.96253311)
        self.x0[ix_['PHI8']] = float(2.0349445303)
        self.x0[ix_['LAM8']] = float(-4.050000443)
        self.x0[ix_['PHI9']] = float(-28.45607804)
        self.x0[ix_['LAM9']] = float(22.430117198)
        self.x0[ix_['PHI10']] = float(16.034035489)
        self.x0[ix_['LAM10']] = float(-17.28050167)
        self.x0[ix_['PHI11']] = float(0.8717052037)
        self.x0[ix_['LAM11']] = float(-0.833052840)
        v_['PHIS'] = v_['DPHI/S']
        v_['START'] = v_['PHI0']+v_['PHIS']
        pass
        for I in range(int(v_['2']),int(v_['S-1'])+1):
            v_['RI'] = float(I)
            v_['PHIS'] = v_['DPHI/S']*v_['RI']
            v_['START'] = v_['PHI0']+v_['PHIS']
            pass
            v_['LAMS'] = v_['DLAM/S']*v_['RI']
            v_['START'] = v_['LAM0']+v_['LAMS']
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE', iet_)
        elftv = loaset(elftv,it,0,'PHI1')
        elftv = loaset(elftv,it,1,'PHI2')
        elftv = loaset(elftv,it,2,'LAM1')
        elftv = loaset(elftv,it,3,'LAM2')
        [it,iet_,_] = s2mpj_ii( 'eE3', iet_)
        elftv = loaset(elftv,it,0,'PHI1')
        elftv = loaset(elftv,it,1,'PHI2')
        elftv = loaset(elftv,it,2,'LAM1')
        elftp = []
        elftp = loaset(elftp,it,0,'LAMF')
        [it,iet_,_] = s2mpj_ii( 'eE2', iet_)
        elftv = loaset(elftv,it,0,'PHI1')
        elftv = loaset(elftv,it,1,'LAM1')
        elftp = loaset(elftp,it,0,'PHIF')
        elftp = loaset(elftp,it,1,'LAMF')
        [it,iet_,_] = s2mpj_ii( 'eE1', iet_)
        elftv = loaset(elftv,it,0,'PHI1')
        elftp = loaset(elftp,it,0,'PHIF')
        elftp = loaset(elftp,it,1,'LAMF')
        elftp = loaset(elftp,it,2,'LAMS')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        ename = 'E0,1'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE1')
        ielftype = arrset(ielftype,ie,iet_["eE1"])
        vname = 'PHI1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='PHIF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PHI0']))
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM0']))
        posep = np.where(elftp[ielftype[ie]]=='LAMS')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM1']))
        ename = 'E0,2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE2')
        ielftype = arrset(ielftype,ie,iet_["eE2"])
        vname = 'PHI2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='PHIF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PHI0']))
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM0']))
        ename = 'E1,2'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE3')
        ielftype = arrset(ielftype,ie,iet_["eE3"])
        vname = 'PHI1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'PHI2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM1']))
        ename = 'E1,3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE3')
        ielftype = arrset(ielftype,ie,iet_["eE3"])
        vname = 'PHI1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'PHI3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM1']))
        ename = 'E2,3'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE')
        ielftype = arrset(ielftype,ie,iet_["eE"])
        vname = 'PHI2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'PHI3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM2')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['4']),int(v_['S-1'])+1):
            v_['I2'] = -2+I
            ename = 'E'+str(int(v_['I2']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE')
            ielftype = arrset(ielftype,ie,iet_["eE"])
            ename = 'E'+str(int(v_['I2']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'PHI'+str(int(v_['I2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I2']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'PHI'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='PHI2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I2']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'LAM'+str(int(v_['I2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I2']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'LAM'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='LAM2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['I1'] = -1+I
            ename = 'E'+str(int(v_['I1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE')
            ielftype = arrset(ielftype,ie,iet_["eE"])
            ename = 'E'+str(int(v_['I1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'PHI'+str(int(v_['I1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'PHI'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='PHI2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'LAM'+str(int(v_['I1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'E'+str(int(v_['I1']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'LAM'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
            posev = np.where(elftv[ielftype[ie]]=='LAM2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'E10,12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE2')
        ielftype = arrset(ielftype,ie,iet_["eE2"])
        vname = 'PHI10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='PHIF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PHI12']))
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM12']))
        ename = 'E11,12'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eE2')
        ielftype = arrset(ielftype,ie,iet_["eE2"])
        vname = 'PHI11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='PHI1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'LAM11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-1000.0),float(1000.0),None)
        posev = np.where(elftv[ielftype[ie]]=='LAM1')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        posep = np.where(elftp[ielftype[ie]]=='PHIF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['PHI12']))
        posep = np.where(elftp[ielftype[ie]]=='LAMF')[0]
        self.elpar = loaset(self.elpar,ie,posep[0],float(v_['LAM12']))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['R0,1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['E0,1'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['2']),int(v_['S'])+1):
            v_['I2'] = -2+I
            ig = ig_['R'+str(int(v_['I2']))+','+str(I)]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['I2']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
            v_['I1'] = -1+I
            ig = ig_['R'+str(int(v_['I1']))+','+str(I)]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['E'+str(int(v_['I1']))+','+str(I)]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SANTA               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-AN-21-23"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S1 = np.sin(EV_[0,0])
        S2 = np.sin(EV_[1,0])
        C1 = np.cos(EV_[0,0])
        C2 = np.cos(EV_[1,0])
        C = np.cos(EV_[2,0]-EV_[3,0])
        S = np.sin(EV_[2,0]-EV_[3,0])
        C1C2S = C1*C2*S
        C1C2C = C1*C2*C
        C1S2S = C1*S2*S
        S1C2S = S1*C2*S
        f_   = S1*S2+C1*C2*C
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1*S2-S1*C2*C
            g_[1] = S1*C2-C1*S2*C
            g_[2] = -C1C2S
            g_[3] = C1C2S
            if nargout>2:
                H_ = np.zeros((4,4))
                H_[0,0] = -S1*S2-C1*C2*C
                H_[1,0] = C1*C2+S1*S2*C
                H_[0,1] = H_[1,0]
                H_[1,1] = -S1*S2-C1*C2*C
                H_[2,0] = S1C2S
                H_[0,2] = H_[2,0]
                H_[2,1] = C1S2S
                H_[1,2] = H_[2,1]
                H_[2,2] = -C1C2C
                H_[0,3] = -S1C2S
                H_[3,0] = H_[0,3]
                H_[1,3] = -C1S2S
                H_[3,1] = H_[1,3]
                H_[2,3] = C1C2C
                H_[3,2] = H_[2,3]
                H_[3,3] = -C1C2C
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S1 = np.sin(EV_[0,0])
        S2 = np.sin(EV_[1,0])
        C1 = np.cos(EV_[0,0])
        C2 = np.cos(EV_[1,0])
        C = np.cos(EV_[2,0]-self.elpar[iel_][0])
        S = np.sin(EV_[2,0]-self.elpar[iel_][0])
        f_   = S1*S2+C1*C2*C
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1*S2-S1*C2*C
            g_[1] = S1*C2-C1*S2*C
            g_[2] = -C1*C2*S
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,0] = -S1*S2-C1*C2*C
                H_[1,0] = C1*C2+S1*S2*C
                H_[0,1] = H_[1,0]
                H_[1,1] = -S1*S2-C1*C2*C
                H_[2,0] = S1*C2*S
                H_[0,2] = H_[2,0]
                H_[2,1] = C1*S2*S
                H_[1,2] = H_[2,1]
                H_[2,2] = -C1*C2*C
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S1 = np.sin(EV_[0,0])
        SF = np.sin(self.elpar[iel_][0])
        C1 = np.cos(EV_[0,0])
        CF = np.cos(self.elpar[iel_][0])
        C = np.cos(EV_[1,0]-self.elpar[iel_][1])
        S = np.sin(EV_[1,0]-self.elpar[iel_][1])
        f_   = S1*SF+C1*CF*C
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1*SF-S1*CF*C
            g_[1] = -C1*CF*S
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,0] = -S1*SF-C1*CF*C
                H_[1,0] = S1*CF*S
                H_[0,1] = H_[1,0]
                H_[1,1] = -C1*CF*C
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S1 = np.sin(EV_[0,0])
        SF = np.sin(self.elpar[iel_][0])
        C1 = np.cos(EV_[0,0])
        CF = np.cos(self.elpar[iel_][0])
        C = np.cos(self.elpar[iel_][2]-self.elpar[iel_][1])
        S = np.sin(self.elpar[iel_][2]-self.elpar[iel_][1])
        f_   = S1*SF+C1*CF*C
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = C1*SF-S1*CF*C
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -S1*SF-C1*CF*C
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

