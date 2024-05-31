from s2xlib import *
class  TAX13322(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TAX13322
#    --------
# 
#    An optimal income tax model with multidimensional taxpayer types,
#    due to Judd, Ma, Saunders & Su
# 
#    Source:
#    Kenneth L. Judd, Ma,  Michael A. Saunders and Che-Lin Su
#    "Optimal Income Taxation with Multidimensional Taxpayer Types"
#    Working Paper, Hoover Institute, Stanford University, 2017
# 
#    SIF input: Nick Gould, July 2018 based on the AMPL model pTAX5Dncl
# 
#    "If ever there was an example that exhibited the stupidity of SIF,
#     this is it. NIMG"
# 
#    classification = "OOR2-MN-72-1261"
# 
#    parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'TAX13322'

    def __init__(self, *args): 
        import numpy as np
        pbm      = structtype()
        pb       = structtype()
        pb.name  = self.name
        pb.sifpbname = 'TAX13322'
        pbm.name = self.name
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['NA'] = int(1);  #  SIF file default value
        else:
            v_['NA'] = int(args[0])
        if nargin<2:
            v_['NB'] = int(3);  #  SIF file default value
        else:
            v_['NB'] = int(args[1])
        if nargin<3:
            v_['NC'] = int(3);  #  SIF file default value
        else:
            v_['NC'] = int(args[2])
        if nargin<4:
            v_['ND'] = int(2);  #  SIF file default value
        else:
            v_['ND'] = int(args[3])
        if nargin<5:
            v_['NE'] = int(2);  #  SIF file default value
        else:
            v_['NE'] = int(args[4])
        v_['0'] = 0
        v_['1'] = 1
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['ONE'] = 1.0e0
        v_['TWO'] = 2.0e0
        v_['THREE'] = 3.0e0
        v_['NBD'] = v_['NB']*v_['ND']
        v_['NCE'] = v_['NC']*v_['NE']
        v_['NP'] = v_['NBD']*v_['NCE']
        v_['NP'] = v_['NP']*v_['NA']
        v_['NPM1'] = -1+v_['NP']
        v_['M'] = v_['NP']*v_['NPM1']
        v_['OMEGA1'] = v_['ONE']/v_['TWO']
        v_['OMEGA2'] = v_['TWO']/v_['THREE']
        v_['THETA1'] = v_['ONE']/v_['THREE']
        v_['THETA2'] = v_['ONE']/v_['TWO']
        v_['THETA3'] = v_['TWO']/v_['THREE']
        v_['PSI1'] = 1.0e0
        v_['PSI2'] = 1.5e0
        v_['W1'] = 2.0e0
        v_['W2'] = 2.5e0
        v_['W3'] = 3.0e0
        v_['W4'] = 3.5e0
        v_['W5'] = 4.0e0
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['LAM'+str(I)+','+str(P)+','+str(Q)] = 1.0e+0
        v_['Q'] = v_['1']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA1']
        v_['Q'] = v_['2']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA2']
        v_['Q'] = v_['3']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA1']
        v_['Q'] = v_['4']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA2']
        v_['Q'] = v_['5']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA1']
        v_['Q'] = v_['6']
        v_['RA'+str(int(v_['Q']))] = v_['ONE']/v_['OMEGA2']
        for I in range(int(v_['1']),int(v_['NA'])+1):
            v_['LOGW'] = np.log(v_['W'+str(I)])
            v_['P'] = v_['1']
            v_['-THETA'] = -1.0*v_['THETA1']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI1']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
            v_['P'] = v_['2']
            v_['-THETA'] = -1.0*v_['THETA1']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI2']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
            v_['P'] = v_['3']
            v_['-THETA'] = -1.0*v_['THETA2']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI1']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
            v_['P'] = v_['4']
            v_['-THETA'] = -1.0*v_['THETA2']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI2']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
            v_['P'] = v_['5']
            v_['-THETA'] = -1.0*v_['THETA3']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI1']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
            v_['P'] = v_['6']
            v_['-THETA'] = -1.0*v_['THETA3']
            v_['RB'] = v_['LOGW']*v_['-THETA']
            v_['RE'] = np.exp(v_['RB'])
            v_['RB'] = v_['RE']*v_['PSI2']
            v_['RB'+str(I)+','+str(int(v_['P']))] = v_['RB']/v_['-THETA']
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = np.array([])
        xscale    = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    [iv,ix_,_] = s2x_ii('C'+str(I)+','+str(P)+','+str(Q),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'C'+str(I)+','+str(P)+','+str(Q))
                    [iv,ix_,_] = s2x_ii('Y'+str(I)+','+str(P)+','+str(Q),ix_)
                    pb.xnames=arrset(pb.xnames,iv,'Y'+str(I)+','+str(P)+','+str(Q))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A       = lil_matrix((1000000,1000000))
        pbm.gscale  = np.array([])
        pbm.grnames = np.array([])
        cnames      = np.array([])
        pb.cnames   = np.array([])
        gtype       = np.array([])
        [ig,ig_,_] = s2x_ii('OBJ',ig_)
        gtype = arrset(gtype,ig,'<>')
        for L in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2x_ii('I'+str(L),ig_)
            gtype = arrset(gtype,ig,'>=')
            cnames = arrset(cnames,ig,'I'+str(L))
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['LAMBDA'] = v_['LAM'+str(I)+','+str(P)+','+str(Q)]
                    v_['-LAMBDA'] = -1.0e0*v_['LAMBDA']
                    [ig,ig_,_] = s2x_ii('T',ig_)
                    gtype = arrset(gtype,ig,'>=')
                    cnames = arrset(cnames,ig,'T')
                    iv = ix_['Y'+str(I)+','+str(P)+','+str(Q)]
                    pbm.A[ig,iv] = float(v_['LAMBDA'])+pbm.A[ig,iv]
                    iv = ix_['C'+str(I)+','+str(P)+','+str(Q)]
                    pbm.A[ig,iv] = float(v_['-LAMBDA'])+pbm.A[ig,iv]
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
        pb.xlower = np.zeros((pb.n,1))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = np.full((pb.n,1),-float('Inf'))
        pb.xupper = np.full((pb.n,1),+float('Inf'))
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    pb.xlower[ix_['C'+str(I)+','+str(P)+','+str(Q)]] = 0.1e0
                    pb.xlower[ix_['Y'+str(I)+','+str(P)+','+str(Q)]] = 0.1e0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = np.full((pb.n,1),float(0.1e0))
        pb.y0 = np.full((pb.m,1),float(0.1e0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2x_ii( 'eA1', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eA2', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eA3', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eA4', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eA5', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eA6', iet_)
        elftv = loaset(elftv,it,0,'C')
        [it,iet_,_] = s2x_ii( 'eB1', iet_)
        elftv = loaset(elftv,it,0,'Y')
        [it,iet_,_] = s2x_ii( 'eB2', iet_)
        elftv = loaset(elftv,it,0,'Y')
        [it,iet_,_] = s2x_ii( 'eB3', iet_)
        elftv = loaset(elftv,it,0,'Y')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        pbm.elftype = np.array([])
        ielftype    = np.array([])
        pbm.elvar   = []
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A1-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA1')
                    ielftype = arrset(ielftype, ie, iet_["eA1"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A2-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA2')
                    ielftype = arrset(ielftype, ie, iet_["eA2"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A3-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA3')
                    ielftype = arrset(ielftype, ie, iet_["eA3"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A4-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA4')
                    ielftype = arrset(ielftype, ie, iet_["eA4"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A5-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA5')
                    ielftype = arrset(ielftype, ie, iet_["eA5"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'A6-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eA6')
                    ielftype = arrset(ielftype, ie, iet_["eA6"])
                    vname = 'C'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='C')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'B1-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eB1')
                    ielftype = arrset(ielftype, ie, iet_["eB1"])
                    vname = 'Y'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'B2-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eB2')
                    ielftype = arrset(ielftype, ie, iet_["eB2"])
                    vname = 'Y'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                for Q in range(int(v_['1']),int(v_['NCE'])+1):
                    ename = 'B3-'+str(I)+','+str(P)+','+str(Q)
                    [ie,ie_,_] = s2x_ii(ename,ie_)
                    pbm.elftype = arrset(pbm.elftype,ie,'eB3')
                    ielftype = arrset(ielftype, ie, iet_["eB3"])
                    vname = 'Y'+str(I)+','+str(P)+','+str(Q)
                    [iv,ix_,pb] = s2x_nlx(vname,ix_,pb,1,None,None,0.1e0)
                    posev = find(elftv[ielftype[ie]],lambda x:x=='Y')
                    pbm.elvar = loaset(pbm.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt   = []
        for ig in np.arange(0,ngrp):
            pbm.grelt.append(np.array([]))
        pbm.grftype = np.array([])
        pbm.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NA'])+1):
            for P in range(int(v_['1']),int(v_['NBD'])+1):
                v_['Q'] = v_['1']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['Q'] = v_['2']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['Q'] = v_['3']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['Q'] = v_['4']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['Q'] = v_['5']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['Q'] = v_['6']
                v_['R'] = (v_['RA'+str(int(v_['Q']))]*v_['LAM'+str(I)+','+str(P)+','+
                     str(int(v_['Q']))])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(P)+','+str(int(v_['Q']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
        for I in range(int(v_['1']),int(v_['NA'])+1):
            v_['LOGW'] = np.log(v_['W'+str(I)])
            for Q in range(int(v_['1']),int(v_['NCE'])+1):
                v_['P'] = v_['1']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                ig = ig_['OBJ']
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['P'] = v_['2']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['P'] = v_['3']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['P'] = v_['4']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['P'] = v_['5']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
                v_['P'] = v_['6']
                v_['R'] = (v_['RB'+str(I)+','+str(int(v_['P']))]*v_['LAM'+str(I)+','+
                     str(int(v_['P']))+','+str(Q)])
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(Q)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['R']))
        v_['L'] = v_['0']
        for I in range(int(v_['1']),int(v_['NA'])+1):
            v_['I+1'] = 1+I
            v_['I-1'] = -1+I
            v_['LOGW'] = np.log(v_['W'+str(I)])
            v_['P'] = v_['1']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['P'] = v_['2']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['P'] = v_['3']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['P'] = v_['4']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['P'] = v_['5']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['P'] = v_['6']
            v_['P+1'] = 1+v_['P']
            v_['P-1'] = -1+v_['P']
            v_['RB'] = v_['RB'+str(I)+','+str(int(v_['P']))]
            v_['-RB'] = -1.0*v_['RB']
            v_['RA'] = v_['RA'+str(int(v_['1']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['P']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A1-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['2']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A1-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['1']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['2']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A2-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['1'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['3']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A2-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['2']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['3']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['2'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['4']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['3']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['4']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A4-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['3'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['5']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A4-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['4']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['5']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A5-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['4'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['6']),int(v_['NCE'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A5-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['5']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            v_['RA'] = v_['RA'+str(int(v_['6']))]
            v_['-RA'] = -1.0*v_['RA']
            for S in range(int(v_['1']),int(v_['P-1'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for S in range(int(v_['P+1']),int(v_['NBD'])+1):
                for R in range(int(v_['1']),int(v_['NA'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['1']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['2']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['3']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['4']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['5']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
                    v_['L'] = v_['L']+v_['1']
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(S)+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['1']),int(v_['I-1'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for R in range(int(v_['I+1']),int(v_['NA'])+1):
                for T in range(int(v_['1']),int(v_['NCE'])+1):
                    v_['L'] = v_['L']+v_['1']
                    ig = ig_['I'+str(int(v_['L']))]
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['A6-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                    posel = len(pbm.grelt[ig])
                    pbm.grelt  = (
                          loaset(pbm.grelt,ig,posel,ie_['B3-'+str(R)+','+str(int(v_['NBD']))+','+str(T)]))
                    nlc = np.union1d(nlc,np.array([ig]))
                    pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
            for T in range(int(v_['1']),int(v_['5'])+1):
                v_['L'] = v_['L']+v_['1']
                ig = ig_['I'+str(int(v_['L']))]
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['A6-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RA']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(int(v_['6']))]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['RB']))
                posel = len(pbm.grelt[ig])
                pbm.grelt  = (
                      loaset(pbm.grelt,ig,posel,ie_['B3-'+str(I)+','+str(int(v_['NBD']))+','+str(T)]))
                nlc = np.union1d(nlc,np.array([ig]))
                pbm.grelw = loaset(pbm.grelw,ig,posel,float(v_['-RB']))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = np.zeros((ngrp,1))
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = np.full((pb.m,1),-float('Inf'))
        pb.cupper = np.full((pb.m,1),+float('Inf'))
        pb.clower[np.arange(pb.nle+pb.neq,pb.m)] = np.zeros((pb.nge,1))
        #%%%%%%%%%%%%%%%%%  RESIZE A %%%%%%%%%%%%%%%%%%%%%%
        pbm.A.resize(ngrp,pb.n)
        pbm.A      = pbm.A.tocsr()
        sA1,sA2    = pbm.A.shape
        pbm.Ashape = [ sA1, sA2 ]
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        lincons =  find(pbm.congrps,lambda x:x in np.setdiff1d(nlc,pbm.congrps))
        pb.pbclass = "OOR2-MN-72-1261"
        self.pb = pb; self.pbm = pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def e_globs(pbm):

        import numpy as np
        pbm.efpar = np.array([]);
        pbm.efpar = arrset( pbm.efpar,0,0.1e0)
        return pbm

    @staticmethod
    def eA1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 0.0e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 0.0e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 1.0e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA4(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 1.0e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA5(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 1.5e0
        OMEGA = 1.0e0/2.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eA6(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        ALPHA = 1.5e0
        OMEGA = 2.0e0/3.0e0
        OM1 = OMEGA-1.0e0
        OM2 = OMEGA-2.0e0
        CMA = EV_[0]-ALPHA
        BIG = CMA>=pbm.efpar[0]
        if BIG!=0:
            F = CMA**OMEGA
        if BIG!=0:
            G = OMEGA*CMA**OM1
        if BIG!=0:
            H = OMEGA*OM1*CMA**OM2
        if BIG==0:
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar[0]**OM1
        if BIG==0:
            C2 = OMEGA*OM1*pbm.efpar[0]**OM2
        if BIG==0:
            F = ((1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar[0]**OMEGA+
             C1*CMA+0.5e0*C2*CMA*CMA)
        if BIG==0:
            G = C1+C2*CMA
        if BIG==0:
            H = C2
        f_   = F
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = G
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = H
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eB1(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        THETA = 1.0e0/3.0e0
        f_   = EV_[0]**THETA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = THETA*EV_[0]**(THETA-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = THETA*(THETA-1.0e0)*EV_[0]**(THETA-2.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eB2(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        THETA = 1.0e0/2.0e0
        f_   = EV_[0]**THETA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = THETA*EV_[0]**(THETA-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = THETA*(THETA-1.0e0)*EV_[0]**(THETA-2.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eB3(pbm,nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        THETA = 2.0e0/3.0e0
        f_   = EV_[0]**THETA
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = THETA*EV_[0]**(THETA-1.0e0)
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = THETA*(THETA-1.0e0)*EV_[0]**(THETA-2.0e0)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

