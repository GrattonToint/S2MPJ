from s2mpjlib import *
class  ROBOT(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This program solves the displacement optimization problem in 
#    REDUNDANT ROBOTS. A redundant robot is one which has more links than 
#    the dimensions it moves in.  Because this redundancy allows almost
#    infinite combinations of joint angles for a particular orientation of
#    end-effector of a robot, choosing an optimum combination has always been a
#    problem of research. 
#    The ROBOT considered here is a 7 link robot moving in 2 dimensional space.
# 
#    Source: an exercize for L. Watson course on LANCELOT in the Spring 1993.
#    B.Benhabib, R.G.Fenton and A.A.Goldberg, 
#    "Analytical trajectory optimization of seven degrees of freedom redundant
#    robot",  
#    Transactions of the Canadian Society for Mechanical Engineering,
#    vol.11(4), 1987, pp 197-200.
# 
#    SIF input: Manish Sabu at Virginia Tech., Spring 1993.
#               Minor modifications by Ph. L. Toint, April 1993.
# 
#    classification = "C-CQOR2-MY-14-2"
# 
#  This segment describes the initial values of angles (by THnIN)
#   and final position of the end effector (by XPOS and YPOS)
#   these values can be changed here according to the needs of the user.
#  The segment also defines the upper and lower bounds of the various joint 
#   angles (by HIGH and DOWN)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ROBOT'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['TH1IN'] = 0.0
        v_['TH2IN'] = 0.0
        v_['TH3IN'] = 0.0
        v_['TH4IN'] = 0.0
        v_['TH5IN'] = 0.0
        v_['TH6IN'] = 0.0
        v_['TH7IN'] = 0.0
        v_['XPOS'] = 4.0
        v_['YPOS'] = 4.0
        v_['HIGH'] = 2.356194
        v_['DOWN'] = -2.356194
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('TH1',ix_)
        self.xnames=arrset(self.xnames,iv,'TH1')
        [iv,ix_,_] = s2mpj_ii('TH2',ix_)
        self.xnames=arrset(self.xnames,iv,'TH2')
        [iv,ix_,_] = s2mpj_ii('TH3',ix_)
        self.xnames=arrset(self.xnames,iv,'TH3')
        [iv,ix_,_] = s2mpj_ii('TH4',ix_)
        self.xnames=arrset(self.xnames,iv,'TH4')
        [iv,ix_,_] = s2mpj_ii('TH5',ix_)
        self.xnames=arrset(self.xnames,iv,'TH5')
        [iv,ix_,_] = s2mpj_ii('TH6',ix_)
        self.xnames=arrset(self.xnames,iv,'TH6')
        [iv,ix_,_] = s2mpj_ii('TH7',ix_)
        self.xnames=arrset(self.xnames,iv,'TH7')
        [iv,ix_,_] = s2mpj_ii('TH1I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH1I')
        [iv,ix_,_] = s2mpj_ii('TH2I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH2I')
        [iv,ix_,_] = s2mpj_ii('TH3I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH3I')
        [iv,ix_,_] = s2mpj_ii('TH4I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH4I')
        [iv,ix_,_] = s2mpj_ii('TH5I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH5I')
        [iv,ix_,_] = s2mpj_ii('TH6I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH6I')
        [iv,ix_,_] = s2mpj_ii('TH7I',ix_)
        self.xnames=arrset(self.xnames,iv,'TH7I')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        [ig,ig_,_] = s2mpj_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        [ig,ig_,_] = s2mpj_ii('CONSTR1',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONSTR1')
        [ig,ig_,_] = s2mpj_ii('CONSTR2',ig_)
        gtype = arrset(gtype,ig,'==')
        cnames = arrset(cnames,ig,'CONSTR2')
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
        self.gconst = arrset(self.gconst,ig_['CONSTR1'],float(v_['XPOS']))
        self.gconst = arrset(self.gconst,ig_['CONSTR2'],float(v_['YPOS']))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        self.xlower[ix_['TH1']] = v_['DOWN']
        self.xupper[ix_['TH1']] = v_['HIGH']
        self.xlower[ix_['TH2']] = v_['DOWN']
        self.xupper[ix_['TH2']] = v_['HIGH']
        self.xlower[ix_['TH3']] = v_['DOWN']
        self.xupper[ix_['TH3']] = v_['HIGH']
        self.xlower[ix_['TH4']] = v_['DOWN']
        self.xupper[ix_['TH4']] = v_['HIGH']
        self.xlower[ix_['TH5']] = v_['DOWN']
        self.xupper[ix_['TH5']] = v_['HIGH']
        self.xlower[ix_['TH6']] = v_['DOWN']
        self.xupper[ix_['TH6']] = v_['HIGH']
        self.xlower[ix_['TH7']] = v_['DOWN']
        self.xupper[ix_['TH7']] = v_['HIGH']
        self.xlower[ix_['TH1I']] = v_['TH1IN']
        self.xupper[ix_['TH1I']] = v_['TH1IN']
        self.xlower[ix_['TH2I']] = v_['TH2IN']
        self.xupper[ix_['TH2I']] = v_['TH2IN']
        self.xlower[ix_['TH3I']] = v_['TH3IN']
        self.xupper[ix_['TH3I']] = v_['TH3IN']
        self.xlower[ix_['TH4I']] = v_['TH4IN']
        self.xupper[ix_['TH4I']] = v_['TH4IN']
        self.xlower[ix_['TH5I']] = v_['TH5IN']
        self.xupper[ix_['TH5I']] = v_['TH5IN']
        self.xlower[ix_['TH6I']] = v_['TH6IN']
        self.xupper[ix_['TH6I']] = v_['TH6IN']
        self.xlower[ix_['TH7I']] = v_['TH7IN']
        self.xupper[ix_['TH7I']] = v_['TH7IN']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['TH1']] = float(0.0)
        self.x0[ix_['TH2']] = float(0.0)
        self.x0[ix_['TH3']] = float(0.0)
        self.x0[ix_['TH4']] = float(0.0)
        self.x0[ix_['TH5']] = float(0.0)
        self.x0[ix_['TH6']] = float(0.0)
        self.x0[ix_['TH7']] = float(0.0)
        self.x0[ix_['TH1I']] = float(0.0)
        self.x0[ix_['TH2I']] = float(0.0)
        self.x0[ix_['TH3I']] = float(0.0)
        self.x0[ix_['TH4I']] = float(0.0)
        self.x0[ix_['TH5I']] = float(0.0)
        self.x0[ix_['TH6I']] = float(0.0)
        self.x0[ix_['TH7I']] = float(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V')
        elftv = loaset(elftv,it,1,'W')
        [it,iet_,_] = s2mpj_ii( 'eCOSTH', iet_)
        elftv = loaset(elftv,it,0,'THETAC')
        [it,iet_,_] = s2mpj_ii( 'eSINTH', iet_)
        elftv = loaset(elftv,it,0,'THETAS')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'TH1SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH1I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH2SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH2I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH3SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH3I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH4SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH4I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH5SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH5I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH6SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH6I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'TH7SQ'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eISQ')
        ielftype = arrset(ielftype,ie,iet_["eISQ"])
        vname = 'TH7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='V')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = 'TH7I'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='W')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C1TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C2TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C3TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C4TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C5TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C6TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'C7TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eCOSTH')
        ielftype = arrset(ielftype,ie,iet_["eCOSTH"])
        vname = 'TH7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAC')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S1TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S2TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S3TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S4TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S5TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S6TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'S7TH'
        [ie,ie_,_] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eSINTH')
        ielftype = arrset(ielftype,ie,iet_["eSINTH"])
        vname = 'TH7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
        posev = np.where(elftv[ielftype[ie]]=='THETAS')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['OBJ']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH1SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH2SQ'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH3SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH4SQ'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH5SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH6SQ'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['TH7SQ'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        ig = ig_['CONSTR1']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C1TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C2TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C3TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C4TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C5TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['C6TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['C7TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.5))
        ig = ig_['CONSTR2']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S1TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['S2TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S3TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['S4TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S5TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,1.)
        posel = posel+1
        self.grelt = loaset(self.grelt,ig,posel,ie_['S6TH'])
        self.grelw = loaset(self.grelw,ig,posel, 1.)
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['S7TH'])
        nlc = np.union1d(nlc,np.array([ig]))
        self.grelw = loaset(self.grelw,ig,posel,float(0.5))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            5.46283877
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CQOR2-MY-14-2"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(self, nargout,*args):

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
    def eCOSTH(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = np.cos(EV_[0,0])
        f_   = TMP
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -np.sin(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -TMP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eSINTH(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        TMP = np.sin(EV_[0,0])
        f_   = TMP
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = np.cos(EV_[0,0])
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -TMP
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

