from s2mpjlib import *
class  BQPGASIM(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BQPGASIM
#    *********
# 
#    The first 50 variable subproblem from BQPGAUSS.
# 
#    SIF input: N. Gould, July 1990.
# 
#    classification = "C-CQBR2-AN-50-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'BQPGASIM'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        [ig,ig_,_] = s2mpj_ii('LINGROUP',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        ngrp   = len(ig_)
        [iv,ix_,_] = s2mpj_ii('1',ix_)
        self.xnames=arrset(self.xnames,iv,'1')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(5.6987e-02))
        [iv,ix_,_] = s2mpj_ii('2',ix_)
        self.xnames=arrset(self.xnames,iv,'2')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-6.1847e-03))
        [iv,ix_,_] = s2mpj_ii('3',ix_)
        self.xnames=arrset(self.xnames,iv,'3')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(5.2516e-03))
        [iv,ix_,_] = s2mpj_ii('4',ix_)
        self.xnames=arrset(self.xnames,iv,'4')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.1729e-02))
        [iv,ix_,_] = s2mpj_ii('5',ix_)
        self.xnames=arrset(self.xnames,iv,'5')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(4.9596e-03))
        [iv,ix_,_] = s2mpj_ii('6',ix_)
        self.xnames=arrset(self.xnames,iv,'6')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.9271e-03))
        [iv,ix_,_] = s2mpj_ii('7',ix_)
        self.xnames=arrset(self.xnames,iv,'7')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.2185e-02))
        [iv,ix_,_] = s2mpj_ii('8',ix_)
        self.xnames=arrset(self.xnames,iv,'8')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.3238e-02))
        [iv,ix_,_] = s2mpj_ii('9',ix_)
        self.xnames=arrset(self.xnames,iv,'9')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.5134e-02))
        [iv,ix_,_] = s2mpj_ii('10',ix_)
        self.xnames=arrset(self.xnames,iv,'10')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.2247e-02))
        [iv,ix_,_] = s2mpj_ii('11',ix_)
        self.xnames=arrset(self.xnames,iv,'11')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(2.3741e-02))
        [iv,ix_,_] = s2mpj_ii('12',ix_)
        self.xnames=arrset(self.xnames,iv,'12')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-9.7666e-02))
        [iv,ix_,_] = s2mpj_ii('13',ix_)
        self.xnames=arrset(self.xnames,iv,'13')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(9.8702e-02))
        [iv,ix_,_] = s2mpj_ii('14',ix_)
        self.xnames=arrset(self.xnames,iv,'14')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(7.8901e-04))
        [iv,ix_,_] = s2mpj_ii('15',ix_)
        self.xnames=arrset(self.xnames,iv,'15')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(5.1663e-04))
        [iv,ix_,_] = s2mpj_ii('16',ix_)
        self.xnames=arrset(self.xnames,iv,'16')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.7477e-04))
        [iv,ix_,_] = s2mpj_ii('17',ix_)
        self.xnames=arrset(self.xnames,iv,'17')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.1795e-03))
        [iv,ix_,_] = s2mpj_ii('18',ix_)
        self.xnames=arrset(self.xnames,iv,'18')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.7351e-02))
        [iv,ix_,_] = s2mpj_ii('19',ix_)
        self.xnames=arrset(self.xnames,iv,'19')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.3439e-03))
        [iv,ix_,_] = s2mpj_ii('20',ix_)
        self.xnames=arrset(self.xnames,iv,'20')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-5.6977e-02))
        [iv,ix_,_] = s2mpj_ii('21',ix_)
        self.xnames=arrset(self.xnames,iv,'21')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.0040e-02))
        [iv,ix_,_] = s2mpj_ii('22',ix_)
        self.xnames=arrset(self.xnames,iv,'22')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-8.3380e-02))
        [iv,ix_,_] = s2mpj_ii('23',ix_)
        self.xnames=arrset(self.xnames,iv,'23')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-3.7526e-03))
        [iv,ix_,_] = s2mpj_ii('24',ix_)
        self.xnames=arrset(self.xnames,iv,'24')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-9.4555e-04))
        [iv,ix_,_] = s2mpj_ii('25',ix_)
        self.xnames=arrset(self.xnames,iv,'25')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.9258e-03))
        [iv,ix_,_] = s2mpj_ii('26',ix_)
        self.xnames=arrset(self.xnames,iv,'26')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.3959e-03))
        [iv,ix_,_] = s2mpj_ii('27',ix_)
        self.xnames=arrset(self.xnames,iv,'27')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.3749e-03))
        [iv,ix_,_] = s2mpj_ii('28',ix_)
        self.xnames=arrset(self.xnames,iv,'28')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.3677e-03))
        [iv,ix_,_] = s2mpj_ii('29',ix_)
        self.xnames=arrset(self.xnames,iv,'29')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-2.7985e-02))
        [iv,ix_,_] = s2mpj_ii('30',ix_)
        self.xnames=arrset(self.xnames,iv,'30')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.8839e-03))
        [iv,ix_,_] = s2mpj_ii('31',ix_)
        self.xnames=arrset(self.xnames,iv,'31')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.2340e-03))
        [iv,ix_,_] = s2mpj_ii('32',ix_)
        self.xnames=arrset(self.xnames,iv,'32')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-6.8139e-04))
        [iv,ix_,_] = s2mpj_ii('33',ix_)
        self.xnames=arrset(self.xnames,iv,'33')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-3.5838e-02))
        [iv,ix_,_] = s2mpj_ii('34',ix_)
        self.xnames=arrset(self.xnames,iv,'34')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-3.4857e-02))
        [iv,ix_,_] = s2mpj_ii('35',ix_)
        self.xnames=arrset(self.xnames,iv,'35')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(2.8724e-03))
        [iv,ix_,_] = s2mpj_ii('36',ix_)
        self.xnames=arrset(self.xnames,iv,'36')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.6625e-02))
        [iv,ix_,_] = s2mpj_ii('37',ix_)
        self.xnames=arrset(self.xnames,iv,'37')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(1.3571e-02))
        [iv,ix_,_] = s2mpj_ii('38',ix_)
        self.xnames=arrset(self.xnames,iv,'38')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-7.2447e-03))
        [iv,ix_,_] = s2mpj_ii('39',ix_)
        self.xnames=arrset(self.xnames,iv,'39')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.6034e-04))
        [iv,ix_,_] = s2mpj_ii('40',ix_)
        self.xnames=arrset(self.xnames,iv,'40')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.6225e-02))
        [iv,ix_,_] = s2mpj_ii('41',ix_)
        self.xnames=arrset(self.xnames,iv,'41')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(2.2034e-05))
        [iv,ix_,_] = s2mpj_ii('42',ix_)
        self.xnames=arrset(self.xnames,iv,'42')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(5.8844e-02))
        [iv,ix_,_] = s2mpj_ii('43',ix_)
        self.xnames=arrset(self.xnames,iv,'43')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(3.0725e-03))
        [iv,ix_,_] = s2mpj_ii('44',ix_)
        self.xnames=arrset(self.xnames,iv,'44')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(2.8227e-03))
        [iv,ix_,_] = s2mpj_ii('45',ix_)
        self.xnames=arrset(self.xnames,iv,'45')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-2.0681e-02))
        [iv,ix_,_] = s2mpj_ii('46',ix_)
        self.xnames=arrset(self.xnames,iv,'46')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-5.4952e-03))
        [iv,ix_,_] = s2mpj_ii('47',ix_)
        self.xnames=arrset(self.xnames,iv,'47')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(6.2552e-04))
        [iv,ix_,_] = s2mpj_ii('48',ix_)
        self.xnames=arrset(self.xnames,iv,'48')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(3.3782e-02))
        [iv,ix_,_] = s2mpj_ii('49',ix_)
        self.xnames=arrset(self.xnames,iv,'49')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-4.8584e-03))
        [iv,ix_,_] = s2mpj_ii('50',ix_)
        self.xnames=arrset(self.xnames,iv,'50')
        icA  = np.append(icA,[iv])
        irA  = np.append(irA,[ig_['LINGROUP']])
        valA = np.append(valA,float(-1.4371e-03))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-0.1)
        self.xupper = np.full((self.n,1),0.1)
        self.xlower[ix_['1']] = -5.4966e-05
        self.xupper[ix_['1']] = 9.9945e-02
        self.xlower[ix_['2']] = -3.9206e-03
        self.xupper[ix_['3']] = 9.9999e-02
        self.xlower[ix_['4']] = -1.0001e-01
        self.xupper[ix_['4']] = 9.9990e-02
        self.xupper[ix_['5']] = 9.9997e-02
        self.xlower[ix_['6']] = -9.9994e-02
        self.xupper[ix_['6']] = 6.1561e-06
        self.xlower[ix_['7']] = -3.9119e-03
        self.xupper[ix_['7']] = 9.9986e-02
        self.xlower[ix_['8']] = -1.0001e-01
        self.xupper[ix_['8']] = 2.5683e-02
        self.xlower[ix_['9']] = -9.9987e-02
        self.xupper[ix_['9']] = 1.0001e-01
        self.xlower[ix_['10']] = -9.9988e-02
        self.xupper[ix_['10']] = 1.0001e-01
        self.xlower[ix_['11']] = -1.0001e-01
        self.xupper[ix_['11']] = 2.8998e-03
        self.xlower[ix_['12']] = -9.9952e-02
        self.xupper[ix_['12']] = 4.7652e-05
        self.xlower[ix_['13']] = -4.5551e-05
        self.xupper[ix_['13']] = 9.9954e-02
        self.xlower[ix_['14']] = -9.9999e-02
        self.xlower[ix_['16']] = -7.2801e-02
        self.xlower[ix_['18']] = -9.9992e-02
        self.xupper[ix_['18']] = 8.3681e-06
        self.xlower[ix_['20']] = -9.9956e-02
        self.xupper[ix_['20']] = 4.3809e-05
        self.xlower[ix_['22']] = -9.9961e-02
        self.xupper[ix_['22']] = 3.9248e-05
        self.xlower[ix_['25']] = -4.1110e-03
        self.xlower[ix_['29']] = -9.6988e-02
        self.xupper[ix_['29']] = 1.0002e-01
        self.xlower[ix_['32']] = -5.8439e-02
        self.xlower[ix_['33']] = -4.5616e-06
        self.xupper[ix_['33']] = 9.9995e-02
        self.xlower[ix_['34']] = -9.9999e-02
        self.xupper[ix_['34']] = 7.3117e-07
        self.xlower[ix_['35']] = -9.9991e-02
        self.xupper[ix_['35']] = 9.3168e-06
        self.xlower[ix_['36']] = -9.9977e-02
        self.xupper[ix_['36']] = 1.0002e-01
        self.xlower[ix_['37']] = -9.9984e-02
        self.xupper[ix_['37']] = 1.5812e-05
        self.xlower[ix_['39']] = -3.9611e-06
        self.xupper[ix_['39']] = 9.9996e-02
        self.xlower[ix_['40']] = -8.8262e-06
        self.xupper[ix_['40']] = 9.9991e-02
        self.xlower[ix_['41']] = -1.0001e-01
        self.xupper[ix_['41']] = 9.9986e-02
        self.xlower[ix_['43']] = -1.9873e-06
        self.xupper[ix_['43']] = 9.9998e-02
        self.xlower[ix_['45']] = -9.9993e-02
        self.xupper[ix_['45']] = 7.4220e-06
        self.xlower[ix_['46']] = -9.9999e-02
        self.xupper[ix_['46']] = 8.2308e-07
        self.xlower[ix_['47']] = -3.0424e-06
        self.xupper[ix_['47']] = 9.9997e-02
        self.xlower[ix_['48']] = -9.9985e-02
        self.xupper[ix_['48']] = 1.5119e-05
        self.xlower[ix_['49']] = -1.0004e-01
        self.xupper[ix_['49']] = 2.4305e-02
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eOFFDIAG', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        [it,iet_,_] = s2mpj_ii( 'eDIAG', iet_)
        elftv = loaset(elftv,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        ename = 'D   1   1'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        self.x0 = np.zeros((self.n,1))
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  11  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  12  12'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '12'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  20  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '21'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '21'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  21  21'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '21'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  29  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  36  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '37'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  36  37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '37'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  37  37'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '37'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  41  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  42'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '42'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  42'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '42'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  42  42'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '42'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   1  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '1'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  49  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   2   2'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  13  13'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '13'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '22'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '22'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  22  22'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '22'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '30'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '30'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  30  30'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '30'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '38'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  36  38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '38'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  38  38'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '38'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   2  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '2'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   3   3'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '23'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '23'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  23  23'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '23'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  43'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '43'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  43'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '43'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  43  43'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '43'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   3  50'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '3'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '50'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  49  50'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '50'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  50  50'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '50'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   4   4'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  14  14'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '14'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '31'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '31'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  31  31'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '31'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  44'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '44'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  44'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '44'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  44  44'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '44'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   4  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '4'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   5   5'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  15  15'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '15'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '24'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '24'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  24  24'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '24'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  32  32'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '32'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '39'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  36  39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '39'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  39  39'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '39'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  45'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '45'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  45'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '45'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  45  45'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '45'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   5  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '5'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   6   6'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  16  16'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '16'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '25'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '25'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  25  25'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '25'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '33'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '33'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  33  33'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '33'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  46'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '46'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  46'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '46'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  46  46'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '46'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   6  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '6'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   7   7'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  17  17'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '17'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '34'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '34'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  34  34'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '34'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   7  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '7'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   8   8'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '40'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  36  40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '40'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  40  40'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '40'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   8  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '8'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D   9   9'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '26'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '26'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  26  26'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '26'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '35'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  29  35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '35'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  35  35'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '35'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O   9  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '9'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  10  10'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  11'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  20'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  29'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '29'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  36'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '36'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  41'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  47'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '47'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  47'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '47'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  47  47'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '47'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  10  49'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '10'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '49'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '18'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  18  18'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '18'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '27'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  27  27'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '27'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  11  19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '11'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '19'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  19  19'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '19'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  20  28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '20'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '28'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  28  28'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '28'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'O  41  48'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        if newelt:
            self.elftype = arrset(self.elftype,ie,'eOFFDIAG')
            ielftype = arrset(ielftype,ie,iet_['eOFFDIAG'])
        vname = '41'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        vname = '48'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='Y')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        ename = 'D  48  48'
        [ie,ie_,newelt] = s2mpj_ii(ename,ie_)
        self.elftype = arrset(self.elftype,ie,'eDIAG')
        ielftype = arrset(ielftype,ie,iet_["eDIAG"])
        vname = '48'
        [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float(-0.1),float(0.1),None)
        posev = np.where(elftv[ielftype[ie]]=='X')[0]
        self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        ig = ig_['LINGROUP']
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   1   1'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9819e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  11  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  12'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  12'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  12  12'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  20  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  21'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  21'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  21  21'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0362e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  29  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.5103e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  36  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  37'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.5140e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  36  37'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  37  37'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  41  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  42'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  42'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  42  42'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   1  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.7537e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  49  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.8331e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   2   2'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9213e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  13'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  13'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9608e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  13  13'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9698e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  22'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  22'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9608e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  22  22'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(8.9945e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  30'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  30'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9608e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  30  30'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.4885e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  38'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.5140e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  36  38'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9608e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  38  38'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5197e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   2  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.7167e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   3   3'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(8.1209e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(8.1463e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  23'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  23'])
        self.grelw = loaset(self.grelw,ig,posel,float(-8.1463e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  23  23'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-7.3536e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-5.3119e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-6.1506e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  43'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  43'])
        self.grelw = loaset(self.grelw,ig,posel,float(-8.1463e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  43  43'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.9480e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   3  50'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.7566e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  49  50'])
        self.grelw = loaset(self.grelw,ig,posel,float(-8.1463e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  50  50'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   4   4'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.8141e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  14'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  14'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.8225e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  14  14'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.8228e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.5487e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  31'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  31'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.8225e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  31  31'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.8370e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.1312e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  44'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  44'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.8225e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  44  44'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   4  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.7539e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   5   5'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.6350e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  15'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  15'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  15  15'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  24'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  24'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  24  24'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.3863e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  32'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  32'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  32  32'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.7205e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  39'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.5140e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  36  39'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  39  39'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.9971e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  45'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  45'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6427e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  45  45'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   5  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.5757e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   6   6'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  16'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  16'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  16  16'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  25'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  25'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  25  25'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.0289e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  33'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  33'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  33  33'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-6.5144e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-7.5509e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  46'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  46'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  46  46'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   6  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.7565e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   7   7'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9320e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  17'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9709e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  17'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9610e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  17  17'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.9631e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(8.9946e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  34'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  34'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.9610e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  34  34'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.4890e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5199e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   7  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.7188e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   8   8'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.7157e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.7417e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-8.7973e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-6.3446e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  40'])
        self.grelw = loaset(self.grelw,ig,posel,float(6.5140e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  36  40'])
        self.grelw = loaset(self.grelw,ig,posel,float(-9.7431e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  40  40'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-7.3586e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   8  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.5052e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D   9   9'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.9055))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.9605))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  26'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  26'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.9604))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  26  26'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.6517))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  35'])
        self.grelw = loaset(self.grelw,ig,posel,float(9.0300e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  29  35'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.9604))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  35  35'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.9168))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.2464))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O   9  49'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.9243))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  10  10'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0624e+03))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  11'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.9135e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  20'])
        self.grelw = loaset(self.grelw,ig,posel,float(2.9241e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  29'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.6379e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  36'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.9046e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  41'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.2065e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  10  47'])
        self.grelw = loaset(self.grelw,ig,posel,float(7.5507e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  47'])
        self.grelw = loaset(self.grelw,ig,posel,float(-2.9232e+01))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  47  47'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  18'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  18  18'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  27'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  27  27'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  11  19'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  19  19'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  20  28'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  28  28'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['O  41  48'])
        self.grelw = loaset(self.grelw,ig,posel,float(-1.0000e+02))
        posel = len(self.grelt[ig])
        self.grelt = loaset(self.grelt,ig,posel,ie_['D  48  48'])
        self.grelw = loaset(self.grelw,ig,posel,float(1.0000e+02))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO BQPGASIM             -5.519814D-5
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-CQBR2-AN-50-0"
        self.x0        = np.zeros((self.n,1))
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eDIAG(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.5*EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eOFFDIAG(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]
            g_[1] = EV_[0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 1.0
                H_[1,0] = H_[0,1]
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

