from s2mpjlib import *
class  VIBRBEAMNE(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A nonlinear least-squares problem arising from laser-Doppler
#    measurements of a vibrating beam.  The data correspond to a simulated
#    experiment where two laser-Doppler velocimeters take measurements
#    at random points along the centreline of the beam.  These measurements
#    consist of a position (x), an incident angle (p) and the magnitude
#    of the velocity along the line of sight (v).
#    The problem is then to fit
# 
#                          2      3                    2     3
#        v = (c + c x + c x  + c x ) cos[ d + d x + d x + d x  - p ]
#              0   1     2      3          0   1     2     3
#            <---- magnitude ----->       <------ phase ----->
# 
#    in the least-squares sense.
# 
#    Source:
#    a modification of an exercize for L. Watson course on LANCELOT in
#    the Spring 1993. Compared to the original proposal, the unnecessary
#    elements were removed as well as an unnecessary constraint on the phase.
# 
#    SIF input: Ph. L. Toint, May 1993, based on a proposal by
#               D. E. Montgomery, Virginia Tech., April 1993.
#    Nonlinear-equations version of VIBRBEAM.SIF, Nick Gould, Jan 2020.
# 
#    classification = "C-CNOR2-MN-8-30"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 7 II 2026
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'VIBRBEAMNE'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['0'] = 0
        v_['1'] = 1
        v_['3'] = 3
        v_['m'] = 30
        v_['x1'] = 39.1722
        v_['x2'] = 53.9707
        v_['x3'] = 47.9829
        v_['x4'] = 12.5925
        v_['x5'] = 16.5414
        v_['x6'] = 18.9548
        v_['x7'] = 27.7168
        v_['x8'] = 31.9201
        v_['x9'] = 45.6830
        v_['x10'] = 22.2524
        v_['x11'] = 33.9805
        v_['x12'] = 6.8425
        v_['x13'] = 35.1677
        v_['x14'] = 33.5682
        v_['x15'] = 43.3659
        v_['x16'] = 13.3835
        v_['x17'] = 25.7273
        v_['x18'] = 21.0230
        v_['x19'] = 10.9755
        v_['x20'] = 1.5323
        v_['x21'] = 45.4416
        v_['x22'] = 14.5431
        v_['x23'] = 22.4313
        v_['x24'] = 29.0144
        v_['x25'] = 25.2675
        v_['x26'] = 15.5095
        v_['x27'] = 9.6297
        v_['x28'] = 8.3009
        v_['x29'] = 30.8694
        v_['x30'] = 43.3299
        v_['v1'] = -1.2026
        v_['v2'] = 1.7053
        v_['v3'] = 0.5410
        v_['v4'] = 1.1477
        v_['v5'] = 1.2447
        v_['v6'] = 0.9428
        v_['v7'] = -0.1360
        v_['v8'] = -0.7542
        v_['v9'] = -0.3396
        v_['v10'] = 0.7057
        v_['v11'] = -0.8509
        v_['v12'] = -0.1201
        v_['v13'] = -1.2193
        v_['v14'] = -1.0448
        v_['v15'] = -0.7723
        v_['v16'] = 0.4342
        v_['v17'] = 0.1154
        v_['v18'] = 0.2868
        v_['v19'] = 0.3558
        v_['v20'] = -0.5090
        v_['v21'] = -0.0842
        v_['v22'] = 0.6021
        v_['v23'] = 0.1197
        v_['v24'] = -0.1827
        v_['v25'] = 0.1806
        v_['v26'] = 0.5395
        v_['v27'] = 0.2072
        v_['v28'] = 0.1466
        v_['v29'] = -0.2672
        v_['v30'] = -0.3038
        v_['p1'] = 2.5736
        v_['p2'] = 2.7078
        v_['p3'] = 2.6613
        v_['p4'] = 2.0374
        v_['p5'] = 2.1553
        v_['p6'] = 2.2195
        v_['p7'] = 2.4077
        v_['p8'] = 2.4772
        v_['p9'] = 2.6409
        v_['p10'] = 2.2981
        v_['p11'] = 2.5073
        v_['p12'] = 1.8380
        v_['p13'] = 2.5236
        v_['p14'] = 2.5015
        v_['p15'] = 2.6186
        v_['p16'] = 0.4947
        v_['p17'] = 0.6062
        v_['p18'] = 0.5588
        v_['p19'] = 0.4772
        v_['p20'] = 0.4184
        v_['p21'] = 0.9051
        v_['p22'] = 0.5035
        v_['p23'] = 0.5723
        v_['p24'] = 0.6437
        v_['p25'] = 0.6013
        v_['p26'] = 0.5111
        v_['p27'] = 0.4679
        v_['p28'] = 0.4590
        v_['p29'] = 0.6666
        v_['p30'] = 0.8630
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [iv,ix_,_] = s2mpj_ii('c0',ix_)
        self.xnames=arrset(self.xnames,iv,'c0')
        [iv,ix_,_] = s2mpj_ii('c1',ix_)
        self.xnames=arrset(self.xnames,iv,'c1')
        [iv,ix_,_] = s2mpj_ii('c2',ix_)
        self.xnames=arrset(self.xnames,iv,'c2')
        [iv,ix_,_] = s2mpj_ii('c3',ix_)
        self.xnames=arrset(self.xnames,iv,'c3')
        [iv,ix_,_] = s2mpj_ii('d0',ix_)
        self.xnames=arrset(self.xnames,iv,'d0')
        [iv,ix_,_] = s2mpj_ii('d1',ix_)
        self.xnames=arrset(self.xnames,iv,'d1')
        [iv,ix_,_] = s2mpj_ii('d2',ix_)
        self.xnames=arrset(self.xnames,iv,'d2')
        [iv,ix_,_] = s2mpj_ii('d3',ix_)
        self.xnames=arrset(self.xnames,iv,'d3')
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for i in range(int(v_['1']),int(v_['m'])+1):
            [ig,ig_,_] = s2mpj_ii('f'+str(i),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'f'+str(i))
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
        for i in range(int(v_['1']),int(v_['m'])+1):
            self.gconst = arrset(self.gconst,ig_['f'+str(i)],float(v_['v'+str(i)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.y0 = np.zeros((self.m,1))
        self.x0[ix_['c0']] = float(-3.5)
        self.x0[ix_['c1']] = float(1.0)
        self.x0[ix_['d0']] = float(1.7)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'efun', iet_)
        elftv = loaset(elftv,it,0,'a0')
        elftv = loaset(elftv,it,1,'a1')
        elftv = loaset(elftv,it,2,'a2')
        elftv = loaset(elftv,it,3,'a3')
        elftv = loaset(elftv,it,4,'b')
        elftp = []
        elftp = loaset(elftp,it,0,'y')
        elftp = loaset(elftp,it,1,'q')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for i in range(int(v_['1']),int(v_['m'])+1):
            for j in range(int(v_['0']),int(v_['3'])+1):
                ename = 'fu'+str(i)+','+str(j)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'efun')
                ielftype = arrset(ielftype,ie,iet_["efun"])
                vname = 'd0'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='a0')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'd1'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='a1')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'd2'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='a2')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'd3'
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='a3')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'c'+str(j)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
                posev = np.where(elftv[ielftype[ie]]=='b')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                posep = np.where(elftp[ielftype[ie]]=='y')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['x'+str(i)]))
                posep = np.where(elftp[ielftype[ie]]=='q')[0]
                self.elpar = loaset(self.elpar,ie,posep[0],float(v_['p'+str(i)]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for i in range(int(v_['1']),int(v_['m'])+1):
            v_['y'] = 1.0
            for j in range(int(v_['0']),int(v_['3'])+1):
                ig = ig_['f'+str(i)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['fu'+str(i)+','+str(j)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['y']))
                v_['y'] = v_['y']*v_['x'+str(i)]
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION             0.15644607137
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CNOR2-MN-8-30"
        self.objderlvl = 2
        self.conderlvl = [2]

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def efun(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        y2 = self.elpar[iel_][0]*self.elpar[iel_][0]
        y3 = self.elpar[iel_][0]*y2
        y4 = y2*y2
        y5 = y2*y3
        y6 = y3*y3
        phi  = (
              EV_[0,0]+self.elpar[iel_][0]*(EV_[1,0]+self.elpar[iel_][0]*(EV_[2,0]+self.elpar[iel_][0]*EV_[3,0]))-self.elpar[iel_][1])
        cosphi = np.cos(phi)
        sinphi = np.sin(phi)
        bcos = EV_[4,0]*cosphi
        bsin = EV_[4,0]*sinphi
        f_   = bcos
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = -bsin
            g_[1] = -bsin*self.elpar[iel_][0]
            g_[2] = -bsin*y2
            g_[3] = -bsin*y3
            g_[4] = cosphi
            if nargout>2:
                H_ = np.zeros((5,5))
                H_[0,0] = -bcos
                H_[0,1] = -bcos*self.elpar[iel_][0]
                H_[1,0] = H_[0,1]
                H_[0,2] = -bcos*y2
                H_[2,0] = H_[0,2]
                H_[0,3] = -bcos*y3
                H_[3,0] = H_[0,3]
                H_[0,4] = -sinphi
                H_[4,0] = H_[0,4]
                H_[1,1] = -bcos*y2
                H_[1,2] = -bcos*y3
                H_[2,1] = H_[1,2]
                H_[1,3] = -bcos*y4
                H_[3,1] = H_[1,3]
                H_[1,4] = -sinphi*self.elpar[iel_][0]
                H_[4,1] = H_[1,4]
                H_[2,2] = -bcos*y4
                H_[2,3] = -bcos*y5
                H_[3,2] = H_[2,3]
                H_[2,4] = -sinphi*y2
                H_[4,2] = H_[2,4]
                H_[3,3] = -bcos*y6
                H_[3,4] = -sinphi*y3
                H_[4,3] = H_[3,4]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

