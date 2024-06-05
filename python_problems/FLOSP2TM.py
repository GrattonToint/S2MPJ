from s2mpjlib import *
class  FLOSP2TM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FLOSP2TM
#    *********
# 
#    A  two-dimensional base  flow  problem in an inclined enclosure.
# 
#    Temperature constant y = +/- 1 boundary conditions
#    Middling Reynold's number
# 
#    The flow is considered in a square of length 2,  centered on the
#    origin and aligned with the x-y axes. The square is divided into
#    4 n ** 2  sub-squares,  each of  length 1 / n.  The differential
#    equation is replaced by  discrete nonlinear equations at each of 
#    the grid points. 
# 
#    The differential equation relates the vorticity, temperature and
#    a stream function.
#    
#    Source: 
#    J. N. Shadid
#    "Experimental and computational study of the stability
#    of Natural convection flow in an inclined enclosure",
#    Ph. D. Thesis, University of Minnesota, 1989,
#    problem SP2 (pp.128-130), 
# 
#    SIF input: Nick Gould, August 1993.
# 
#    classification = "NQR2-MY-V-V"
# 
#    Half the number of discretization intervals
#    Number of variables = 3(2M+1)**2 
# 
#           Alternative values for the SIF file parameters:
# IE M                   1              $-PARAMETER n=27
# IE M                   2              $-PARAMETER n=75
# IE M                   5              $-PARAMETER n=363     original value
# IE M                   8              $-PARAMETER n=867
# IE M                   10             $-PARAMETER n=1323
# IE M                   15             $-PARAMETER n=2883
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'FLOSP2TM'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'FLOSP2TM'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['M'] = int(1);  #  SIF file default value
        else:
            v_['M'] = int(args[0])
        if nargin<2:
            v_['RA'] = float(1.0e+5);  #  SIF file default value
        else:
            v_['RA'] = float(args[1])
        v_['PI/4'] = np.arctan(1.0)
        v_['PI'] = 4.0*v_['PI/4']
        v_['AX'] = 1.0
        v_['THETA'] = 0.5*v_['PI']
        v_['A1'] = 0.0
        v_['A2'] = 1.0
        v_['A3'] = 0.0
        v_['B1'] = 0.0
        v_['B2'] = 1.0
        v_['B3'] = 1.0
        v_['F1'] = 1.0
        v_['F2'] = 0.0
        v_['F3'] = 0.0
        v_['G1'] = 1.0
        v_['G2'] = 0.0
        v_['G3'] = 0.0
        v_['M-1'] = -1+v_['M']
        v_['-M'] = -1*v_['M']
        v_['-M+1'] = -1*v_['M-1']
        v_['1/H'] = float(v_['M'])
        v_['-1/H'] = -1.0*v_['1/H']
        v_['2/H'] = 2.0*v_['1/H']
        v_['-2/H'] = -2.0*v_['1/H']
        v_['H'] = 1.0/v_['1/H']
        v_['H2'] = v_['H']*v_['H']
        v_['1/H2'] = v_['1/H']*v_['1/H']
        v_['-2/H2'] = -2.0*v_['1/H2']
        v_['1/2H'] = 0.5*v_['1/H']
        v_['-1/2H'] = -0.5*v_['1/H']
        v_['AXX'] = v_['AX']*v_['AX']
        v_['SINTHETA'] = np.sin(v_['THETA'])
        v_['COSTHETA'] = np.cos(v_['THETA'])
        v_['PI1'] = v_['AX']*v_['RA']
        v_['PI1'] = v_['PI1']*v_['COSTHETA']
        v_['PI1'] = -0.5*v_['PI1']
        v_['-PI1'] = -1.0*v_['PI1']
        v_['PI2'] = v_['AXX']*v_['RA']
        v_['PI2'] = v_['PI2']*v_['SINTHETA']
        v_['PI2'] = 0.5*v_['PI2']
        v_['-PI2'] = -1.0*v_['PI2']
        v_['2A1'] = 2.0*v_['A1']
        v_['2B1'] = 2.0*v_['B1']
        v_['2F1'] = 2.0*v_['F1']
        v_['2G1'] = 2.0*v_['G1']
        v_['2F1/AX'] = v_['2F1']/v_['AX']
        v_['2G1/AX'] = v_['2G1']/v_['AX']
        v_['AX/2'] = 0.5*v_['AX']
        v_['AXX/2'] = 0.5*v_['AXX']
        v_['AXX/4'] = 0.25*v_['AXX']
        v_['2AX'] = 2.0*v_['AX']
        v_['2AXX'] = 2.0*v_['AXX']
        v_['2/AX'] = 2.0/v_['AX']
        v_['2/AXH'] = v_['2/H']/v_['AX']
        v_['-2/AXH'] = -1.0*v_['2/AXH']
        v_['PI1/2H'] = v_['PI1']*v_['1/2H']
        v_['-PI1/2H'] = v_['PI1']*v_['-1/2H']
        v_['PI2/2H'] = v_['PI2']*v_['1/2H']
        v_['-PI2/2H'] = v_['PI2']*v_['-1/2H']
        v_['2A1/H'] = v_['2A1']*v_['1/H']
        v_['-2A1/H'] = v_['2A1']*v_['-1/H']
        v_['2B1/H'] = v_['2B1']*v_['1/H']
        v_['-2B1/H'] = v_['2B1']*v_['-1/H']
        v_['2F1/AXH'] = v_['2F1/AX']*v_['1/H']
        v_['-2F1/AXH'] = v_['2F1/AX']*v_['-1/H']
        v_['2G1/AXH'] = v_['2G1/AX']*v_['1/H']
        v_['-2G1/AXH'] = v_['2G1/AX']*v_['-1/H']
        v_['AX/H2'] = v_['AX']*v_['1/H2']
        v_['-AX/H2'] = -1.0*v_['AX/H2']
        v_['AX/4H2'] = 0.25*v_['AX/H2']
        v_['-AX/4H2'] = -0.25*v_['AX/H2']
        v_['AXX/H2'] = v_['AXX']*v_['1/H2']
        v_['-2AXX/H2'] = -2.0*v_['AXX/H2']
        v_['1'] = 1
        v_['2'] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for J in range(int(v_['-M']),int(v_['M'])+1):
            for I in range(int(v_['-M']),int(v_['M'])+1):
                [iv,ix_,_] = s2mpj_ii('OM'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'OM'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('PH'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'PH'+str(I)+','+str(J))
                [iv,ix_,_] = s2mpj_ii('PS'+str(I)+','+str(J),ix_)
                pb.xnames=arrset(pb.xnames,iv,'PS'+str(I)+','+str(J))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        for J in range(int(v_['-M+1']),int(v_['M-1'])+1):
            v_['J+'] = 1+J
            v_['J-'] = -1+J
            for I in range(int(v_['-M+1']),int(v_['M-1'])+1):
                v_['I+'] = 1+I
                v_['I-'] = -1+I
                [ig,ig_,_] = s2mpj_ii('S'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'S'+str(I)+','+str(J))
                iv = ix_['OM'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(int(v_['I+']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(int(v_['I-']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(I)+','+str(int(v_['J+']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(I)+','+str(int(v_['J-']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(int(v_['I+']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['-PI1/2H'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(int(v_['I-']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['PI1/2H'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(I)+','+str(int(v_['J+']))]
                pbm.A[ig,iv] = float(v_['-PI2/2H'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(I)+','+str(int(v_['J-']))]
                pbm.A[ig,iv] = float(v_['PI2/2H'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('V'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'V'+str(I)+','+str(J))
                iv = ix_['PS'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/H2'])+pbm.A[ig,iv]
                iv = ix_['PS'+str(int(v_['I+']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['PS'+str(int(v_['I-']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['PS'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['PS'+str(I)+','+str(int(v_['J+']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['PS'+str(I)+','+str(int(v_['J-']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['OM'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['AXX/4'])+pbm.A[ig,iv]
                [ig,ig_,_] = s2mpj_ii('E'+str(I)+','+str(J),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'E'+str(I)+','+str(J))
                iv = ix_['PH'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(int(v_['I+']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(int(v_['I-']))+','+str(J)]
                pbm.A[ig,iv] = float(v_['1/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(I)+','+str(J)]
                pbm.A[ig,iv] = float(v_['-2AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(I)+','+str(int(v_['J+']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
                iv = ix_['PH'+str(I)+','+str(int(v_['J-']))]
                pbm.A[ig,iv] = float(v_['AXX/H2'])+pbm.A[ig,iv]
        for K in range(int(v_['-M']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['M']))]
            pbm.A[ig,iv] = float(v_['2A1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['M-1']))]
            pbm.A[ig,iv] = float(v_['-2A1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['M']))]
            pbm.A[ig,iv] = float(v_['A2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['-M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['-M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['-M+1']))]
            pbm.A[ig,iv] = float(v_['2B1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['-M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['-M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['-M']))]
            pbm.A[ig,iv] = float(v_['-2B1/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(K)+','+str(int(v_['-M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(K)+','+str(int(v_['-M'])))
            iv = ix_['PH'+str(K)+','+str(int(v_['-M']))]
            pbm.A[ig,iv] = float(v_['B2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['2F1/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['M-1']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['-2F1/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['F2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['-M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['-M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['-M+1']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['2G1/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['-M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['-M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['-M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['-2G1/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('T'+str(int(v_['-M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'T'+str(int(v_['-M']))+','+str(K))
            iv = ix_['PH'+str(int(v_['-M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['G2'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(K)+','+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(K)+','+str(int(v_['M'])))
            iv = ix_['PS'+str(K)+','+str(int(v_['M']))]
            pbm.A[ig,iv] = float(v_['-2/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(K)+','+str(int(v_['M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(K)+','+str(int(v_['M'])))
            iv = ix_['PS'+str(K)+','+str(int(v_['M-1']))]
            pbm.A[ig,iv] = float(v_['2/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(K)+','+str(int(v_['-M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(K)+','+str(int(v_['-M'])))
            iv = ix_['PS'+str(K)+','+str(int(v_['-M+1']))]
            pbm.A[ig,iv] = float(v_['2/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(K)+','+str(int(v_['-M'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(K)+','+str(int(v_['-M'])))
            iv = ix_['PS'+str(K)+','+str(int(v_['-M']))]
            pbm.A[ig,iv] = float(v_['-2/H'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(int(v_['M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(int(v_['M']))+','+str(K))
            iv = ix_['PS'+str(int(v_['M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['-2/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(int(v_['M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(int(v_['M']))+','+str(K))
            iv = ix_['PS'+str(int(v_['M-1']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['2/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(int(v_['-M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(int(v_['-M']))+','+str(K))
            iv = ix_['PS'+str(int(v_['-M+1']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['2/AXH'])+pbm.A[ig,iv]
            [ig,ig_,_] = s2mpj_ii('V'+str(int(v_['-M']))+','+str(K),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'V'+str(int(v_['-M']))+','+str(K))
            iv = ix_['PS'+str(int(v_['-M']))+','+str(K)]
            pbm.A[ig,iv] = float(v_['-2/AXH'])+pbm.A[ig,iv]
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
        for K in range(int(v_['-M']),int(v_['M'])+1):
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['T'+str(K)+','+str(int(v_['M']))],float(v_['A3'])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['T'+str(K)+','+str(int(v_['-M']))],float(v_['B3'])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['T'+str(int(v_['M']))+','+str(K)],float(v_['F3'])))
            pbm.gconst  = (
                  arrset(pbm.gconst,ig_['T'+str(int(v_['-M']))+','+str(K)],float(v_['G3'])))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for K in range(int(v_['-M']),int(v_['M'])+1):
            pb.xlower[ix_['PS'+str(K)+','+str(int(v_['-M']))]] = 1.0
            pb.xupper[ix_['PS'+str(K)+','+str(int(v_['-M']))]] = 1.0
            pb.xlower[ix_['PS'+str(int(v_['-M']))+','+str(K)]] = 1.0
            pb.xupper[ix_['PS'+str(int(v_['-M']))+','+str(K)]] = 1.0
            pb.xlower[ix_['PS'+str(K)+','+str(int(v_['M']))]] = 1.0
            pb.xupper[ix_['PS'+str(K)+','+str(int(v_['M']))]] = 1.0
            pb.xlower[ix_['PS'+str(int(v_['M']))+','+str(K)]] = 1.0
            pb.xupper[ix_['PS'+str(int(v_['M']))+','+str(K)]] = 1.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'ePROD', iet_)
        elftv = loaset(elftv,it,0,'PSIM')
        elftv = loaset(elftv,it,1,'PSIP')
        elftv = loaset(elftv,it,2,'PHIM')
        elftv = loaset(elftv,it,3,'PHIP')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for J in range(int(v_['-M+1']),int(v_['M-1'])+1):
            v_['J+'] = 1+J
            v_['J-'] = -1+J
            for I in range(int(v_['-M+1']),int(v_['M-1'])+1):
                v_['I+'] = 1+I
                v_['I-'] = -1+I
                ename = 'E'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                pb.x0 = np.zeros((pb.n,1))
                vname = 'PS'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PSIP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PS'+str(I)+','+str(int(v_['J-']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PSIM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PH'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PHIP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PH'+str(int(v_['I-']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PHIM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                ename = 'F'+str(I)+','+str(J)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                pbm.elftype = arrset(pbm.elftype,ie,'ePROD')
                ielftype = arrset(ielftype, ie, iet_["ePROD"])
                vname = 'PS'+str(int(v_['I+']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PSIP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PS'+str(int(v_['I-']))+','+str(J)
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PSIM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PH'+str(I)+','+str(int(v_['J+']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PHIP')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
                vname = 'PH'+str(I)+','+str(int(v_['J-']))
                [iv,ix_,pb] = s2mpj_nlx(vname,ix_,pb,1,None,None,None)
                posev = find(elftv[ielftype[ie]],lambda x:x=='PHIM')
                pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for J in range(int(v_['-M+1']),int(v_['M-1'])+1):
            for I in range(int(v_['-M+1']),int(v_['M-1'])+1):
                ig = ig_['E'+str(I)+','+str(J)]
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['E'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-AX/4H2']))
                posel = len(pbm.grelt[ig])
                pbm.grelt = loaset(pbm.grelt,ig,posel,ie_['F'+str(I)+','+str(J)])
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['AX/4H2']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
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
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "NQR2-MY-V-V"
        pb.x0          = np.zeros((pb.n,1))
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def ePROD(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,4))
        IV_ = np.zeros(2)
        U_[0,1] = U_[0,1]+1
        U_[0,0] = U_[0,0]-1
        U_[1,3] = U_[1,3]+1
        U_[1,2] = U_[1,2]-1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]
            g_[1] = IV_[0]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

