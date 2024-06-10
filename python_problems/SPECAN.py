from s2mpjlib import *
class  SPECAN(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SPECAN
#    *********
# 
#    Source: a problem in spectral analysis suggested
#    by J. Eriksson and P. Lindstrom in "A Parallel Algorithm
#    for Bound Constrained Nonlinear Least Squares", UMEA TR S-901 87
# 
#    SIF input: Michael Ferris, July 1993
# 
#    classification = "SBR2-AN-V-0"
# 
#    Number of Gaussians
# 
#           Alternative values for the SIF file parameters:
# IE K                   1              $-PARAMETER
# IE K                   2              $-PARAMETER
# IE K                   3              $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'SPECAN'

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
            v_['K'] = int(3);  #  SIF file default value
        else:
            v_['K'] = int(args[0])
        v_['N'] = 3
        v_['M'] = 5000
        v_['RealM'] = float(v_['M'])
        v_['H'] = 25.0/v_['RealM']
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['ONE'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for p in range(int(v_['1']),int(v_['K'])+1):
            for j in range(int(v_['1']),int(v_['N'])+1):
                [iv,ix_,_] = s2mpj_ii('X'+str(p)+','+str(j),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(p)+','+str(j))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for p in range(int(v_['1']),int(v_['K'])+1):
            for I in range(int(v_['1']),int(v_['M'])+1):
                [ig,ig_,_] = s2mpj_ii('OBJ'+str(p)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        v_['SOLN1,1'] = 19.0
        v_['SOLN1,2'] = 4.2
        v_['SOLN1,3'] = 1.2
        v_['SOLN2,1'] = 8.0
        v_['SOLN2,2'] = 2.5
        v_['SOLN2,3'] = 4.6
        v_['SOLN3,1'] = 10.0
        v_['SOLN3,2'] = 2.0
        v_['SOLN3,3'] = 2.6
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['RI'] = float(I)
            v_['IH'] = v_['H']*v_['RI']
            v_['TI'] = v_['ONE']+v_['IH']
            v_['Differ'] = v_['TI']-v_['SOLN1,2']
            v_['Numer'] = v_['Differ']*v_['Differ']
            v_['Denom'] = v_['SOLN1,3']*v_['SOLN1,3']
            v_['Differ'] = v_['Numer']/v_['Denom']
            v_['Ratio'] = 0.0-v_['Differ']
            v_['ERat'] = np.exp(v_['Ratio'])
            v_['Yi1'] = v_['SOLN1,1']*v_['ERat']
            v_['Differ'] = v_['TI']-v_['SOLN2,2']
            v_['Numer'] = v_['Differ']*v_['Differ']
            v_['Denom'] = v_['SOLN2,3']*v_['SOLN2,3']
            v_['Differ'] = v_['Numer']/v_['Denom']
            v_['Ratio'] = 0.0-v_['Differ']
            v_['ERat'] = np.exp(v_['Ratio'])
            v_['Yi2'] = v_['SOLN2,1']*v_['ERat']
            v_['Differ'] = v_['TI']-v_['SOLN3,2']
            v_['Numer'] = v_['Differ']*v_['Differ']
            v_['Denom'] = v_['SOLN3,3']*v_['SOLN3,3']
            v_['Differ'] = v_['Numer']/v_['Denom']
            v_['Ratio'] = 0.0-v_['Differ']
            v_['ERat'] = np.exp(v_['Ratio'])
            v_['Yi3'] = v_['SOLN3,1']*v_['ERat']
            for p in range(int(v_['1']),int(v_['K'])+1):
                pbm.gconst  = (
                      arrset(pbm.gconst,ig_['OBJ'+str(p)+','+str(I)],float(v_['Yi'+str(p)])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),float('inf'))
        v_['LOWER1,1'] = 15.0
        v_['LOWER1,2'] = 3.5
        v_['LOWER1,3'] = 0.3
        v_['LOWER2,1'] = 5.0
        v_['LOWER2,2'] = 2.2
        v_['LOWER2,3'] = 2.6
        v_['LOWER3,1'] = 5.0
        v_['LOWER3,2'] = 1.2
        v_['LOWER3,3'] = 1.3
        v_['UPPER1,1'] = 31.0
        v_['UPPER1,2'] = 6.3
        v_['UPPER1,3'] = 3.7
        v_['UPPER2,1'] = 15.0
        v_['UPPER2,2'] = 5.3
        v_['UPPER2,3'] = 6.2
        v_['UPPER3,1'] = 14.0
        v_['UPPER3,2'] = 3.3
        v_['UPPER3,3'] = 2.8
        for p in range(int(v_['1']),int(v_['K'])+1):
            for j in range(int(v_['1']),int(v_['N'])+1):
                pb.xlower[ix_['X'+str(p)+','+str(j)]] = v_['LOWER'+str(p)+','+str(j)]
                pb.xupper[ix_['X'+str(p)+','+str(j)]] = v_['UPPER'+str(p)+','+str(j)]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        v_['START1,1'] = 25.0
        v_['START1,2'] = 5.2
        v_['START1,3'] = 3.2
        v_['START2,1'] = 7.0
        v_['START2,2'] = 4.1
        v_['START2,3'] = 3.6
        v_['START3,1'] = 11.6
        v_['START3,2'] = 1.9
        v_['START3,3'] = 2.2
        for p in range(int(v_['1']),int(v_['K'])+1):
            for j in range(int(v_['1']),int(v_['N'])+1):
                pb.x0[ix_['X'+str(p)+','+str(j)]] = float(v_['START'+str(p)+','+str(j)])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eEXPSQ', iet_)
        elftv = loaset(elftv,it,0,'U')
        elftv = loaset(elftv,it,1,'V')
        elftv = loaset(elftv,it,2,'W')
        elftp = []
        elftp = loaset(elftp,it,0,'T')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for p in range(int(v_['1']),int(v_['K'])+1):
            for I in range(int(v_['1']),int(v_['M'])+1):
                ename = 'E'+str(p)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eEXPSQ')
                ielftype = arrset(ielftype, ie, iet_["eEXPSQ"])
                vname = 'X'+str(p)+','+str(int(v_['1']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='U')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(p)+','+str(int(v_['2']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(p)+','+str(int(v_['3']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='W')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                v_['RI'] = float(I)
                v_['IH'] = v_['H']*v_['RI']
                v_['TI'] = v_['ONE']+v_['IH']
                posep = find(elftp[ielftype[ie]],lambda x:x=='T')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],float(v_['TI']))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gSQUARE',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            pbm.grftype = arrset(pbm.grftype,ig,'gSQUARE')
        for p in range(int(v_['1']),int(v_['K'])+1):
            for I in range(int(v_['1']),int(v_['M'])+1):
                ig = ig_['OBJ'+str(p)+','+str(I)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(p)+','+str(I)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "SBR2-AN-V-0"
        self.pb = pb; self.pbm = pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eEXPSQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        R = (pbm.elpar[iel_][0]-EV_[1])**2
        S = EV_[2]**2
        E = np.exp(-R/S)
        f_   = EV_[0]*E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E
            g_[1] = 2.0*(pbm.elpar[iel_][0]-EV_[1])*EV_[0]*E/S
            g_[2] = 2.0*R*EV_[0]*E/(S*EV_[2])
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = 2.0*(pbm.elpar[iel_][0]-EV_[1])*E/S
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0*R*E/(S*EV_[2])
                H_[2,0] = H_[0,2]
                H_[1,1] = (2.0*EV_[0]*E/S)*(2.0*R/S-1.0)
                H_[1,2] = 4.0*(pbm.elpar[iel_][0]-EV_[1])*EV_[0]*E/(S*EV_[2])*(R/S-1.0)
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*R*EV_[0]*E/(S**3)*(2.0*R-3.0*S)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gSQUARE(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= 5.0e-1*GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 1.0e+0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

