from s2xlib import *
class  MINSURFO(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINSURFO
#    *********
# 
#    Find the surface with minimal area, given boundary conditions, 
#    and above an obstacle.
# 
#    This is problem 17 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, December 2000
# 
#    classification = "OBR2-AN-V-V"
# 
#  grid points in x direction (fixed at 50 in COPS)
# 
# IE NX                  25
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MINSURFO'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'MINSURFO'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['NX'] = 50
        v_['NY'] = 100
        v_['0'] = 0
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['NX+1'] = 1+v_['NX']
        v_['NY+1'] = 1+v_['NY']
        v_['RNX+1'] = float(v_['NX+1'])
        v_['RNY+1'] = float(v_['NY+1'])
        v_['HX'] = 1.0/v_['RNX+1']
        v_['HY'] = 1.0/v_['RNY+1']
        v_['AREA'] = v_['HX']*v_['HY']
        v_['AREA'] = 0.5*v_['AREA']
        v_['1/AREA'] = 1.0/v_['AREA']
        v_['1/HX'] = 1.0/v_['HX']
        v_['1/HX2'] = v_['1/HX']*v_['1/HX']
        v_['1/HY'] = 1.0/v_['HY']
        v_['1/HY2'] = v_['1/HY']*v_['1/HY']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['0']),int(v_['NX+1'])+1):
            for J in range(int(v_['0']),int(v_['NY+1'])+1):
                [iv,ix_,_] = s2x_ii('V'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'V'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                [ig,ig_,_] = s2x_ii('A'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/AREA']))
        for I in range(int(v_['1']),int(v_['NX+1'])+1):
            for J in range(int(v_['1']),int(v_['NY+1'])+1):
                [ig,ig_,_] = s2x_ii('B'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<>')
                pbm.gscale = arrset(pbm.gscale,ig,float(v_['1/AREA']))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = len(ix_)
        ngrp   = len(ig_)
        pbm.objgrps = np.arange(ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['A'+str(I)+','+str(J)],float(-1.0))
        for I in range(int(v_['1']),int(v_['NX+1'])+1):
            for J in range(int(v_['1']),int(v_['NY+1'])+1):
                pbm.gconst = arrset(pbm.gconst,ig_['B'+str(I)+','+str(J)],float(-1.0))
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        v_['1/4HX'] = 0.25/v_['HX']
        v_['3/4HX'] = 0.75/v_['HX']
        v_['1/4HY'] = 0.25/v_['HY']
        v_['3/4HY'] = 0.75/v_['HY']
        v_['3/4HX'] = 0.9999999999+v_['3/4HX']
        v_['3/4HY'] = 0.9999999999+v_['3/4HY']
        v_['1/4HX'] = int(np.fix(v_['1/4HX']))
        v_['1/4HY'] = int(np.fix(v_['1/4HY']))
        v_['3/4HX'] = int(np.fix(v_['3/4HX']))
        v_['3/4HY'] = int(np.fix(v_['3/4HY']))
        for I in range(int(v_['1/4HX']),int(v_['3/4HX'])+1):
            for J in range(int(v_['1/4HY']),int(v_['3/4HY'])+1):
                pb.xlower[ix_['V'+str(I)+','+str(J)]] = 1.0
        for J in range(int(v_['0']),int(v_['NY+1'])+1):
            pb.xlower[ix_['V'+str(int(v_['0']))+','+str(J)]] = 0.0
            pb.xupper[ix_['V'+str(int(v_['0']))+','+str(J)]] = 0.0
            pb.xlower[ix_['V'+str(int(v_['NX+1']))+','+str(J)]] = 0.0
            pb.xupper[ix_['V'+str(int(v_['NX+1']))+','+str(J)]] = 0.0
        for I in range(int(v_['0']),int(v_['NX+1'])+1):
            v_["I"] = float(I)
            v_['VIJ'] = 2.0*I
            v_['VIJ'] = v_['VIJ']*v_['HX']
            v_['VIJ'] = -1.0+v_['VIJ']
            v_['VIJ'] = v_['VIJ']*v_['VIJ']
            v_['VIJ'] = v_['ONE']-v_['VIJ']
            pb.xlower[ix_['V'+str(I)+','+str(int(v_['0']))]] = v_['VIJ']
            pb.xupper[ix_['V'+str(I)+','+str(int(v_['0']))]] = v_['VIJ']
            pb.xlower[ix_['V'+str(I)+','+str(int(v_['NY+1']))]] = v_['VIJ']
            pb.xupper[ix_['V'+str(I)+','+str(int(v_['NY+1']))]] = v_['VIJ']
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        for I in range(int(v_['0']),int(v_['NX+1'])+1):
            v_["I"] = float(I)
            v_['VIJ'] = 2.0*I
            v_['VIJ'] = v_['VIJ']*v_['HX']
            v_['VIJ'] = -1.0+v_['VIJ']
            v_['VIJ'] = v_['VIJ']*v_['VIJ']
            v_['VIJ'] = v_['ONE']-v_['VIJ']
            for J in range(int(v_['0']),int(v_['NY+1'])+1):
                pb.x0[ix_['V'+str(I)+','+str(J)]] = float(v_['VIJ'])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eISQ', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['0']),int(v_['NX'])+1):
            v_['I+1'] = 1+I
            for J in range(int(v_['0']),int(v_['NY'])+1):
                v_['J+1'] = 1+J
                ename = 'I'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'V'+str(int(v_['I+1']))+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'J'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
                ielftype = arrset(ielftype, ie, iet_["eISQ"])
                vname = 'V'+str(I)+','+str(int(v_['J+1']))
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'V'+str(I)+','+str(J)
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for J in range(int(v_['0']),int(v_['NY+1'])+1):
            v_['J1'] = 1+J
            ename = 'J'+str(int(v_['NX+1']))+','+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
            ielftype = arrset(ielftype, ie, iet_["eISQ"])
            ename = 'J'+str(int(v_['NX+1']))+','+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'V'+str(int(v_['NX+1']))+','+str(int(v_['J1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'J'+str(int(v_['NX+1']))+','+str(J)
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'V'+str(int(v_['NX+1']))+','+str(J)
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['0']),int(v_['NX+1'])+1):
            v_['I1'] = 1+I
            ename = 'I'+str(I)+','+str(int(v_['NY+1']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'eISQ')
            ielftype = arrset(ielftype, ie, iet_["eISQ"])
            ename = 'I'+str(I)+','+str(int(v_['NY+1']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'V'+str(int(v_['I1']))+','+str(int(v_['NY+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V1')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'I'+str(I)+','+str(int(v_['NY+1']))
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'V'+str(I)+','+str(int(v_['NY+1']))
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='V2')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2x_ii('gROOT',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['0']),int(v_['NX'])+1):
            for J in range(int(v_['0']),int(v_['NY'])+1):
                ig = ig_['A'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gROOT')
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['I'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/HX2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['J'+str(I)+','+str(J)])
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/HY2']))
        for I in range(int(v_['1']),int(v_['NX+1'])+1):
            v_['I-1'] = -1+I
            for J in range(int(v_['1']),int(v_['NY+1'])+1):
                v_['J-1'] = -1+J
                ig = ig_['B'+str(I)+','+str(J)]
                pbm.grftype = arrset(pbm.grftype,ig,'gROOT')
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['I'+str(int(v_['I-1']))+','+str(J)]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/HX2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['J'+str(I)+','+str(int(v_['J-1']))]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['1/HY2']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( pbm, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        pb.pbclass = "OBR2-AN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eISQ(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((1,2))
        IV_ = np.zeros(1)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        f_   = IV_[0]*IV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gROOT(pbm,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        ROOTAL = np.sqrt(GVAR_)
        f_= ROOTAL
        if nargout>1:
            g_ = 0.5/ROOTAL
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = -0.25/ROOTAL**3
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

