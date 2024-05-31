from s2xlib import *
class  ORBIT2(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#        A reformulation of the discretized optimal control problem
#        ORBIT
#        Consider minimising the time needed for a spacecraft to
#        move from one circular orbit around the earth to another.
# 
#        This problem can be put in the following form:
#        let (x1,x2,x3) be the position and (x4,x5,x6) the velocity
#        let (u1,u2,u3) be the driving force vector
#        let q be the time required, then our problem becomes:
#        MINIMISE        q
#        with the equations of motion:
#                        dx1/dt = hv*q*x4
#                        dx2/dt = hv*q*x5
#                        dx3/dt = hv*q*x6
#                        dx4/dt = hv*q*(hg*x1/r^3-hf*u(1)/m)     (1)
#                        dx5/dt = hv*q*(hg*x1/r^3-hf*u(2)/m)             
#                        dx6/dt = hv*q*(hg*x1/r^3-hf*u(3)/m)
#        
#        with    m = m0-hm*q*t   ( - the mass variation)
#                        
#                r = sqrt( x1^2 + x2^2 + x3^2 )  -  dist. from the center
#                                                   of the earth
#        't' is a rescaled time varying between 0 and 1, and
#        'hv,hf,hg,hm' are scaling constants.
#        the driving force is bounded by 
#                        u1^2 + u2^2 + u3^2 <= 1                 (2)
#        (the rather arbitrary no. '1' representing the max. power
#        of the spacecraft).
#        We choose the initial conditions:
#                x1 = x2 = 0   ,   x3 = 1   -  initial position
#                x5 = x6 = 0   ,   x4 = Vorb    -   corresponding orbital
#                                                   speed
#        and the final conditions are:
#        x1^2 + x2^2 + x3^2 = Rf^2  -   final orbit radius
#        x4^2 + x5^2 + x6^2 = Vforb^2 -  corresponding orbital
#                                                 speed
#        x1*x4  + x2*x5 + x3*x6 = 0  -  direction must be parallel
#                                        to the velocity
#        we have chosen the constants hg,hf,hv,hm so that the x,u vectors
#        are of order one. These correspond to an initial orbit at 150 km
#        above the earth's surface, and a final orbit at 250 km above the 
#        earth's surface.
#        The reduction to an NLP is done in that same way as for CAR.SIF
#        The time taken should be:
#                        q = 315 secs
# 
#    SIF input: Thomas Felici, Universite de Nancy, France, 
#               October 1993, with modifications to compute derivatives
#               by Nick Gould
#      SAVEs removed December 3rd 2014
# 
#    classification = "LOR1-RN-V-V"
# 
# *******************************************
# 
#        M = Number of time nodes.
# 
#        Change this for different resolution
# 
# IE M                              3   $-PARAMETER n= 25, m= 18
# IE M                             10   $-PARAMETER n= 88, m= 67
# IE M                             30   $-PARAMETER n=268, m=207 original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'ORBIT2'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'ORBIT2'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(300);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
        v_['NEQ'] = 6
        v_['NCON'] = 1
        v_['NLIM'] = 3
        v_['NLINEQ'] = 0
        v_['NLINCON'] = 0
        v_['NLINLIM'] = 0
        v_['NTEQ'] = v_['NEQ']+v_['NLINEQ']
        v_['NTCON'] = v_['NCON']+v_['NLINCON']
        v_['NTLIM'] = v_['NLIM']+v_['NLINLIM']
        v_['NX'] = 6
        v_['NU'] = 3
        v_['NQ'] = 1
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['M-1'] = v_['M']-v_['1']
        v_['RM'] = v_['M-1']
        for I in range(int(v_['1']),int(v_['M'])+1):
            v_['I-1'] = I-v_['1']
            v_['RI1'] = v_['I-1']
            v_['T'+str(I)] = v_['RI1']/v_['RM']
        v_['RT'] = 6371.0e0
        v_['RT*RT'] = v_['RT']*v_['RT']
        v_['MU'] = 9.81e-3*v_['RT*RT']
        v_['ONE'] = 1.0e+0
        v_['PI/4'] = np.arctan(v_['ONE'])
        v_['PI'] = 4.0e+0*v_['PI/4']
        v_['R0'] = 1.5e+2+v_['RT']
        v_['R0*R0'] = v_['R0']*v_['R0']
        v_['MU/R0'] = v_['MU']/v_['R0']
        v_['VORB'] = np.sqrt(v_['MU/R0'])
        v_['RF'] = 2.5e+2+v_['RT']
        v_['M0'] = 3.0e+3
        v_['QT'] = 3.333e+0
        v_['VG'] = 3.0e+0
        v_['TS'] = 1.0e+0
        v_['HG'] = v_['VORB']*v_['R0*R0']
        v_['HG'] = v_['MU']/v_['HG']
        v_['HG'] = v_['TS']*v_['HG']
        v_['HV'] = v_['VORB']/v_['R0']
        v_['HV'] = v_['TS']*v_['HV']
        v_['HF'] = v_['VORB']*v_['M0']
        v_['HF'] = v_['VG']/v_['HF']
        v_['HF'] = v_['QT']*v_['HF']
        v_['HF'] = v_['TS']*v_['HF']
        v_['HM'] = v_['QT']/v_['M0']
        v_['HM'] = v_['TS']*v_['HM']
        v_['VF'] = v_['MU']/v_['RF']
        v_['VF'] = np.sqrt(v_['VF'])
        v_['VF'] = v_['VF']/v_['VORB']
        v_['VFVF'] = v_['VF']*v_['VF']
        v_['RF'] = v_['RF']/v_['R0']
        v_['RFRF'] = v_['RF']*v_['RF']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['NX'])+1):
                [iv,ix_,_] = s2x_ii('X'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'X'+str(I)+','+str(J))
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            for J in range(int(v_['1']),int(v_['NU'])+1):
                [iv,ix_,_] = s2x_ii('U'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'U'+str(I)+','+str(J))
        for J in range(int(v_['1']),int(v_['NQ'])+1):
            [iv,ix_,_] = s2x_ii('Q'+str(J),ix_)
            pb.xnames=arrset(pb.xnames,iv,'Q'+str(J))
            xscale = arrset(xscale,iv,1.0e+2)
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        iv = ix_['Q'+str(v_['1'])]
        pbm.A[ig,iv] = 1.0+pbm.A[ig,iv]
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            for J in range(int(v_['1']),int(v_['NTEQ'])+1):
                [ig,ig_,_] = s2x_ii('K'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'K'+str(I)+','+str(J))
        for J in range(int(v_['1']),int(v_['NTLIM'])+1):
            [ig,ig_,_] = s2x_ii('L'+str(J),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'L'+str(J))
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['NTCON'])+1):
                [ig,ig_,_] = s2x_ii('G'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'<=')
                cnames = arrset(cnames,ig,'G'+str(I)+','+str(J))
                pbm.gscale = arrset(pbm.gscale,ig,1.0e+2)
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
        pbm.gconst = arrset(pbm.gconst,ig_['L'+str(v_['1'])],v_['RFRF'])
        pbm.gconst = arrset(pbm.gconst,ig_['L'+str(v_['3'])],v_['VFVF'])
        for I in range(int(v_['1']),int(v_['M'])+1):
            pbm.gconst = arrset(pbm.gconst,ig_['G'+str(I)+','+str(v_['1'])],1.0e+0)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('inf'))
        pb.xupper = np.full((pb.n,1),float('inf'))
        pb.xlower[ix_['Q'+str(v_['1'])]] = .10000E+03
        pb.xupper[ix_['Q'+str(v_['1'])]] = .40000E+26
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['1'])]] = .00000E+00
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['1'])]] = .00000E+00
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['2'])]] = .00000E+00
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['2'])]] = .00000E+00
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['3'])]] = .10000E+01
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['3'])]] = .10000E+01
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['4'])]] = .10000E+01
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['4'])]] = .10000E+01
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['5'])]] = .00000E+00
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['5'])]] = .00000E+00
        pb.xlower[ix_['X'+str(v_['1'])+','+str(v_['6'])]] = .00000E+00
        pb.xupper[ix_['X'+str(v_['1'])+','+str(v_['6'])]] = .00000E+00
        for I in range(int(v_['2']),int(v_['M'])+1):
            pb.xlower[ix_['X'+str(I)+','+str(v_['1'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['1'])]] = .10000E+26
            pb.xlower[ix_['X'+str(I)+','+str(v_['2'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['2'])]] = .10000E+26
            pb.xlower[ix_['X'+str(I)+','+str(v_['3'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['3'])]] = .10000E+26
            pb.xlower[ix_['X'+str(I)+','+str(v_['4'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['4'])]] = .10000E+26
            pb.xlower[ix_['X'+str(I)+','+str(v_['5'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['5'])]] = .10000E+26
            pb.xlower[ix_['X'+str(I)+','+str(v_['6'])]] = -.10000E+26
            pb.xupper[ix_['X'+str(I)+','+str(v_['6'])]] = .10000E+26
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            pb.xlower[ix_['U'+str(I)+','+str(v_['1'])]] = -.10000E+01
            pb.xupper[ix_['U'+str(I)+','+str(v_['1'])]] = .10000E+01
            pb.xlower[ix_['U'+str(I)+','+str(v_['2'])]] = -.10000E+01
            pb.xupper[ix_['U'+str(I)+','+str(v_['2'])]] = .10000E+01
            pb.xlower[ix_['U'+str(I)+','+str(v_['3'])]] = -.10000E+01
            pb.xupper[ix_['U'+str(I)+','+str(v_['3'])]] = .10000E+01
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.zeros((pb.n,1))
        pb.y0 = np.zeros((pb.m,1))
        if('Q'+str(v_['1']) in ix_):
            pb.x0[ix_['Q'+str(v_['1'])]] = .10000E+03
        else:
            pb.y0  = (
                  arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['Q'+str(v_['1'])]),.10000E+03))
        for I in range(int(v_['1']),int(v_['M'])+1):
            if('X'+str(I)+','+str(v_['1']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['1'])]] = .00000E+00
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['1'])]),.00000E+00))
            if('X'+str(I)+','+str(v_['2']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['2'])]] = .00000E+00
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['2'])]),.00000E+00))
            if('X'+str(I)+','+str(v_['3']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['3'])]] = .10000E+01
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['3'])]),.10000E+01))
            if('X'+str(I)+','+str(v_['4']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['4'])]] = .10000E+01
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['4'])]),.10000E+01))
            if('X'+str(I)+','+str(v_['5']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['5'])]] = .00000E+00
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['5'])]),.00000E+00))
            if('X'+str(I)+','+str(v_['6']) in ix_):
                pb.x0[ix_['X'+str(I)+','+str(v_['6'])]] = .00000E+00
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['X'+str(I)+','+str(v_['6'])]),.00000E+00))
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            if('U'+str(I)+','+str(v_['1']) in ix_):
                pb.x0[ix_['U'+str(I)+','+str(v_['1'])]] = .10000E+01
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['U'+str(I)+','+str(v_['1'])]),.10000E+01))
            if('U'+str(I)+','+str(v_['2']) in ix_):
                pb.x0[ix_['U'+str(I)+','+str(v_['2'])]] = .10000E+01
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['U'+str(I)+','+str(v_['2'])]),.10000E+01))
            if('U'+str(I)+','+str(v_['3']) in ix_):
                pb.x0[ix_['U'+str(I)+','+str(v_['3'])]] = .10000E+01
            else:
                pb.y0  = (
                      arrset(pb.y0,findfirst(pbm.congrps,lambda x:x==ig_['U'+str(I)+','+str(v_['3'])]),.10000E+01))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'SQR', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2x_ii( 'PROD', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2x_ii( 'Ktyp', iet_)
        elftv = loaset(elftv,it,0,'XI1')
        elftv = loaset(elftv,it,1,'XF1')
        elftv = loaset(elftv,it,2,'XI2')
        elftv = loaset(elftv,it,3,'XF2')
        elftv = loaset(elftv,it,4,'XI3')
        elftv = loaset(elftv,it,5,'XF3')
        elftv = loaset(elftv,it,6,'XI4')
        elftv = loaset(elftv,it,7,'XF4')
        elftv = loaset(elftv,it,8,'XI5')
        elftv = loaset(elftv,it,9,'XF5')
        elftv = loaset(elftv,it,10,'XI6')
        elftv = loaset(elftv,it,11,'XF6')
        elftv = loaset(elftv,it,12,'UV1')
        elftv = loaset(elftv,it,13,'UV2')
        elftv = loaset(elftv,it,14,'UV3')
        elftv = loaset(elftv,it,15,'QV1')
        elftp = []
        elftp = loaset(elftp,it,0,'TN')
        elftp = loaset(elftp,it,1,'TN1')
        elftp = loaset(elftp,it,2,'ID')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        pbm.elpar   = []
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            v_['S'] = 1+I
            for J in range(int(v_['1']),int(v_['NTEQ'])+1):
                v_['ReJ'] = J
                ename = 'Ke'+str(I)+','+str(J)
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'Ktyp')
                ielftype = arrset(ielftype, ie, iet_["Ktyp"])
                posep = find(elftp[ielftype[ie]],lambda x:x=='TN')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['T'+str(I)])
                posep = find(elftp[ielftype[ie]],lambda x:x=='TN1')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['T'+str(v_['S'])])
                posep = find(elftp[ielftype[ie]],lambda x:x=='ID')
                pbm.elpar = loaset(pbm.elpar,ie,posep[0],v_['ReJ'])
                vname = 'X'+str(I)+','+str(v_['1'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['1'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(v_['2'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['2'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(v_['3'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['3'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(v_['4'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['4'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF4')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(v_['5'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI5')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['5'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF5')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(I)+','+str(v_['6'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XI6')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'X'+str(v_['S'])+','+str(v_['6'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='XF6')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(v_['1'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='UV1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(v_['2'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='UV2')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'U'+str(I)+','+str(v_['3'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='UV3')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'Q'+str(v_['1'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='QV1')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['NTCON'])+1):
            for I in range(int(v_['1']),int(v_['M-1'])+1):
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['1'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'SQR')
                ielftype = arrset(ielftype, ie, iet_["SQR"])
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['1'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                vname = 'U'+str(I)+','+str(v_['1'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['2'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'SQR')
                ielftype = arrset(ielftype, ie, iet_["SQR"])
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['2'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                vname = 'U'+str(I)+','+str(v_['2'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['3'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'SQR')
                ielftype = arrset(ielftype, ie, iet_["SQR"])
                ename = 'Ge'+str(I)+','+str(J)+','+str(v_['3'])
                [ie,ie_,_] = s2x_ii(ename,ie_)
                vname = 'U'+str(I)+','+str(v_['3'])
                [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='X')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for J in range(int(v_['1']),int(v_['NTCON'])+1):
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['1'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'SQR')
            ielftype = arrset(ielftype, ie, iet_["SQR"])
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['1'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'U'+str(v_['M-1'])+','+str(v_['1'])
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['2'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'SQR')
            ielftype = arrset(ielftype, ie, iet_["SQR"])
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['2'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'U'+str(v_['M-1'])+','+str(v_['2'])
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['3'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            pbm.elftype = arrset(pbm.elftype,ie,'SQR')
            ielftype = arrset(ielftype, ie, iet_["SQR"])
            ename = 'Ge'+str(v_['M'])+','+str(J)+','+str(v_['3'])
            [ie,ie_,_] = s2x_ii(ename,ie_)
            vname = 'U'+str(v_['M-1'])+','+str(v_['3'])
            [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
            posev = find(elftv[ielftype[ie]],lambda x:x=='X')
            pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['1'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['1'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['1'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['1'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['1'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['2'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['1'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['1'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['3'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PROD')
        ielftype = arrset(ielftype, ie, iet_["PROD"])
        ename = 'Le'+str(v_['2'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['1'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['4'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PROD')
        ielftype = arrset(ielftype, ie, iet_["PROD"])
        ename = 'Le'+str(v_['2'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['2'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['5'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'PROD')
        ielftype = arrset(ielftype, ie, iet_["PROD"])
        ename = 'Le'+str(v_['2'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['3'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['2'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['6'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['3'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['3'])+','+str(v_['1'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['4'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['3'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['3'])+','+str(v_['2'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['5'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        ename = 'Le'+str(v_['3'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        pbm.elftype = arrset(pbm.elftype,ie,'SQR')
        ielftype = arrset(ielftype, ie, iet_["SQR"])
        ename = 'Le'+str(v_['3'])+','+str(v_['3'])
        [ie,ie_,_] = s2x_ii(ename,ie_)
        vname = 'X'+str(v_['M'])+','+str(v_['6'])
        [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,None)
        posev = find(elftv[ielftype[ie]],lambda x:x=='X')
        pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['M-1'])+1):
            for J in range(int(v_['1']),int(v_['NTEQ'])+1):
                ig = ig_['K'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Ke'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['M'])+1):
            for J in range(int(v_['1']),int(v_['NTCON'])+1):
                ig = ig_['G'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['Ge'+str(I)+','+str(J)+','+str(v_['1'])]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['Ge'+str(I)+','+str(J)+','+str(v_['2'])]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['Ge'+str(I)+','+str(J)+','+str(v_['3'])]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NTLIM'])+1):
            ig = ig_['L'+str(I)]
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Le'+str(I)+','+str(v_['1'])])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Le'+str(I)+','+str(v_['2'])])
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
            posel = len(pbm.grelt[ig])
            pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['Le'+str(I)+','+str(v_['3'])])
            nlc = np.union1d(nlc,np.array([ig]))
            pbm.grelw = loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.cupper[np.arange(pb.nle)] = np.zeros((pb.nle,1))
        pb.clower[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        pb.cupper[np.arange(pb.nle,pb.nle+pb.neq)] = np.zeros((pb.neq,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        lxs = len(xscale);
        for j in np.arange(0,min(sA2,pb.n,len(xscale))):
            if not xscale[j] is None and xscale[j] != 0.0 and xscale[j] != 1.0:
                for i in find(pbm.A[:,j],lambda x:x!=0):
                      pbm.A[i,j] = pbm.A[i,j]/xscale[j]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "LOR1-RN-V-V"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
