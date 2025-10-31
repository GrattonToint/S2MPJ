from s2mpjlib import *
class  NCB20(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NCB20
#    *********
# 
#    A banded problem with semi-bandwidth 20.  This problem exhibits frequent
#    negative curvature in the exact Hessian.
# 
#    Source:
#    Ph. Toint, private communication, 1992.
# 
#    SIF input: Ph. Toint, October 1992.
# 
#    classification = "C-COUR2-AN-V-0"
# 
#    Problem dimension
# 
#           Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'NCB20'

    def __init__(self, *args): 
        import numpy as np
        from scipy.sparse import csr_matrix
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        if nargin<1:
            v_['N'] = int(25);  #  SIF file default value
        else:
            v_['N'] = int(args[0])
        v_['P'] = 20
        v_['NY'] = 10
        v_['COND'] = 1.0e4
        v_['1/COND'] = 1.0/v_['COND']
        v_['RP'] = float(v_['P'])
        v_['CL'] = -4.0/v_['RP']
        v_['-P'] = -1*v_['P']
        v_['N-P'] = v_['N']+v_['-P']
        v_['1'] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('X'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'X'+str(I))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            [iv,ix_,_] = s2mpj_ii('Y'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'Y'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        for I in range(int(v_['1']),int(v_['N-P'])+1):
            for J in range(int(v_['1']),int(v_['P'])+1):
                v_['I+J'] = I+J
                v_['I+J-1'] = -1+v_['I+J']
                [ig,ig_,_] = s2mpj_ii('O'+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['X'+str(int(v_['I+J-1']))]])
                valA = np.append(valA,float(v_['CL']))
        [ig,ig_,_] = s2mpj_ii('W',ig_)
        gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        self.gconst = np.full((ngrp,1),-2.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        for I in range(int(v_['1']),int(v_['NY'])+1):
            self.x0[ix_['Y'+str(I)]] = float(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eBP', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftv = loaset(elftv,it,3,'V4')
        elftv = loaset(elftv,it,4,'V5')
        elftv = loaset(elftv,it,5,'V6')
        elftv = loaset(elftv,it,6,'V7')
        elftv = loaset(elftv,it,7,'V8')
        elftv = loaset(elftv,it,8,'V9')
        elftv = loaset(elftv,it,9,'V10')
        elftv = loaset(elftv,it,10,'V11')
        elftv = loaset(elftv,it,11,'V12')
        elftv = loaset(elftv,it,12,'V13')
        elftv = loaset(elftv,it,13,'V14')
        elftv = loaset(elftv,it,14,'V15')
        elftv = loaset(elftv,it,15,'V16')
        elftv = loaset(elftv,it,16,'V17')
        elftv = loaset(elftv,it,17,'V18')
        elftv = loaset(elftv,it,18,'V19')
        elftv = loaset(elftv,it,19,'V20')
        [it,iet_,_] = s2mpj_ii( 'eQR', iet_)
        elftv = loaset(elftv,it,0,'XX')
        [it,iet_,_] = s2mpj_ii( 'en3P', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['N-P'])+1):
            v_['I+1'] = 1+I
            v_['I+2'] = 2+I
            v_['I+3'] = 3+I
            v_['I+4'] = 4+I
            v_['I+5'] = 5+I
            v_['I+6'] = 6+I
            v_['I+7'] = 7+I
            v_['I+8'] = 8+I
            v_['I+9'] = 9+I
            v_['I+10'] = 10+I
            v_['I+11'] = 11+I
            v_['I+12'] = 12+I
            v_['I+13'] = 13+I
            v_['I+14'] = 14+I
            v_['I+15'] = 15+I
            v_['I+16'] = 16+I
            v_['I+17'] = 17+I
            v_['I+18'] = 18+I
            v_['I+19'] = 19+I
            ename = 'E'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eBP')
            ielftype = arrset(ielftype,ie,iet_["eBP"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+1']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+2']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+3']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V4')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+4']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V5')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+5']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V6')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+6']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V7')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+7']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V8')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+8']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V9')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+9']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V10')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+10']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V11')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+11']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V12')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+12']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V13')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+13']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V14')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+14']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V15')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+15']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V16')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+16']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V17')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+17']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V18')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+18']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V19')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'X'+str(int(v_['I+19']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='V20')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['N'])+1):
            ename = 'S'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eQR')
            ielftype = arrset(ielftype,ie,iet_["eQR"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='XX')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NY'])+1):
            ename = 'L'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'en3P')
            ielftype = arrset(ielftype,ie,iet_["en3P"])
            vname = 'X'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            v_['NY+I'] = v_['NY']+I
            vname = 'X'+str(int(v_['NY+I']))
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'Y'+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,float('-Inf'),float('Inf'),None)
            posev = np.where(elftv[ielftype[ie]]=='Z')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['N-P'])+1):
            v_['RI'] = float(I)
            v_['1/I'] = 10.0/v_['RI']
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['E'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['1/I']))
        for I in range(int(v_['1']),int(v_['N'])+1):
            ig = ig_['O'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['S'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        for I in range(int(v_['1']),int(v_['NY'])+1):
            ig = ig_['W']
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['L'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['1/COND']))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass   = "C-COUR2-AN-V-0"
        self.objderlvl = 2


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eBP(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        S = 1.0
        D1 = S+EV_[0]*EV_[0]
        Y1 = EV_[0]/D1
        DY1 = (1.0-2.0*EV_[0]*EV_[0]/D1)/D1
        D2Y1 = (8.0*EV_[0]**3/D1-6.0*EV_[0])/D1**2
        D2 = S+EV_[1]*EV_[1]
        Y2 = EV_[1]/D2
        DY2 = (1.0-2.0*EV_[1]*EV_[1]/D2)/D2
        D2Y2 = (8.0*EV_[1]**3/D2-6.0*EV_[1])/D2**2
        D3 = S+EV_[2]*EV_[2]
        Y3 = EV_[2]/D3
        DY3 = (1.0-2.0*EV_[2]*EV_[2]/D3)/D3
        D2Y3 = (8.0*EV_[2]**3/D3-6.0*EV_[2])/D3**2
        D4 = S+EV_[3]*EV_[3]
        Y4 = EV_[3]/D4
        DY4 = (1.0-2.0*EV_[3]*EV_[3]/D4)/D4
        D2Y4 = (8.0*EV_[3]**3/D4-6.0*EV_[3])/D4**2
        D5 = S+EV_[4]*EV_[4]
        Y5 = EV_[4]/D5
        DY5 = (1.0-2.0*EV_[4]*EV_[4]/D5)/D5
        D2Y5 = (8.0*EV_[4]**3/D5-6.0*EV_[4])/D5**2
        D6 = S+EV_[5]*EV_[5]
        Y6 = EV_[5]/D6
        DY6 = (1.0-2.0*EV_[5]*EV_[5]/D6)/D6
        D2Y6 = (8.0*EV_[5]**3/D6-6.0*EV_[5])/D6**2
        D7 = S+EV_[6]*EV_[6]
        Y7 = EV_[6]/D7
        DY7 = (1.0-2.0*EV_[6]*EV_[6]/D7)/D7
        D2Y7 = (8.0*EV_[6]**3/D7-6.0*EV_[6])/D7**2
        D8 = S+EV_[7]*EV_[7]
        Y8 = EV_[7]/D8
        DY8 = (1.0-2.0*EV_[7]*EV_[7]/D8)/D8
        D2Y8 = (8.0*EV_[7]**3/D8-6.0*EV_[7])/D8**2
        D9 = S+EV_[8]*EV_[8]
        Y9 = EV_[8]/D9
        DY9 = (1.0-2.0*EV_[8]*EV_[8]/D9)/D9
        D2Y9 = (8.0*EV_[8]**3/D9-6.0*EV_[8])/D9**2
        D10 = S+EV_[9]*EV_[9]
        Y10 = EV_[9]/D10
        DY10 = (1.0-2.0*EV_[9]*EV_[9]/D10)/D10
        D2Y10 = (8.0*EV_[9]**3/D10-6.0*EV_[9])/D10**2
        D11 = S+EV_[10]*EV_[10]
        Y11 = EV_[10]/D11
        DY11 = (1.0-2.0*EV_[10]*EV_[10]/D11)/D11
        D2Y11 = (8.0*EV_[10]**3/D11-6.0*EV_[10])/D11**2
        D12 = S+EV_[11]*EV_[11]
        Y12 = EV_[11]/D12
        DY12 = (1.0-2.0*EV_[11]*EV_[11]/D12)/D12
        D2Y12 = (8.0*EV_[11]**3/D12-6.0*EV_[11])/D12**2
        D13 = S+EV_[12]*EV_[12]
        Y13 = EV_[12]/D13
        DY13 = (1.0-2.0*EV_[12]*EV_[12]/D13)/D13
        D2Y13 = (8.0*EV_[12]**3/D13-6.0*EV_[12])/D13**2
        D14 = S+EV_[13]*EV_[13]
        Y14 = EV_[13]/D14
        DY14 = (1.0-2.0*EV_[13]*EV_[13]/D14)/D14
        D2Y14 = (8.0*EV_[13]**3/D14-6.0*EV_[13])/D14**2
        D15 = S+EV_[14]*EV_[14]
        Y15 = EV_[14]/D15
        DY15 = (1.0-2.0*EV_[14]*EV_[14]/D15)/D15
        D2Y15 = (8.0*EV_[14]**3/D15-6.0*EV_[14])/D15**2
        D16 = S+EV_[15]*EV_[15]
        Y16 = EV_[15]/D16
        DY16 = (1.0-2.0*EV_[15]*EV_[15]/D16)/D16
        D2Y16 = (8.0*EV_[15]**3/D16-6.0*EV_[15])/D16**2
        D17 = S+EV_[16]*EV_[16]
        Y17 = EV_[16]/D17
        DY17 = (1.0-2.0*EV_[16]*EV_[16]/D17)/D17
        D2Y17 = (8.0*EV_[16]**3/D17-6.0*EV_[16])/D17**2
        D18 = S+EV_[17]*EV_[17]
        Y18 = EV_[17]/D18
        DY18 = (1.0-2.0*EV_[17]*EV_[17]/D18)/D18
        D2Y18 = (8.0*EV_[17]**3/D18-6.0*EV_[17])/D18**2
        D19 = S+EV_[18]*EV_[18]
        Y19 = EV_[18]/D19
        DY19 = (1.0-2.0*EV_[18]*EV_[18]/D19)/D19
        D2Y19 = (8.0*EV_[18]**3/D19-6.0*EV_[18])/D19**2
        D20 = S+EV_[19]*EV_[19]
        Y20 = EV_[19]/D20
        DY20 = (1.0-2.0*EV_[19]*EV_[19]/D20)/D20
        D2Y20 = (8.0*EV_[19]**3/D20-6.0*EV_[19])/D20**2
        SUM = Y1+Y2+Y3+Y4+Y5+Y19+Y20+Y6+Y7+Y8+Y9+Y10+Y11+Y12+Y13+Y14+Y15+Y16+Y17+Y18
        TWOSUM = SUM+SUM
        f_   = SUM*SUM
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = TWOSUM*DY1
            g_[1] = TWOSUM*DY2
            g_[2] = TWOSUM*DY3
            g_[3] = TWOSUM*DY4
            g_[4] = TWOSUM*DY5
            g_[5] = TWOSUM*DY6
            g_[6] = TWOSUM*DY7
            g_[7] = TWOSUM*DY8
            g_[8] = TWOSUM*DY9
            g_[9] = TWOSUM*DY10
            g_[10] = TWOSUM*DY11
            g_[11] = TWOSUM*DY12
            g_[12] = TWOSUM*DY13
            g_[13] = TWOSUM*DY14
            g_[14] = TWOSUM*DY15
            g_[15] = TWOSUM*DY16
            g_[16] = TWOSUM*DY17
            g_[17] = TWOSUM*DY18
            g_[18] = TWOSUM*DY19
            g_[19] = TWOSUM*DY20
            if nargout>2:
                H_ = np.zeros((20,20))
                H_[0,0] = 2.0*DY1*DY1+TWOSUM*D2Y1
                H_[0,1] = 2.0*DY1*DY2
                H_[1,0] = H_[0,1]
                H_[0,2] = 2.0*DY1*DY3
                H_[2,0] = H_[0,2]
                H_[0,3] = 2.0*DY1*DY4
                H_[3,0] = H_[0,3]
                H_[0,4] = 2.0*DY1*DY5
                H_[4,0] = H_[0,4]
                H_[0,5] = 2.0*DY1*DY6
                H_[5,0] = H_[0,5]
                H_[0,6] = 2.0*DY1*DY7
                H_[6,0] = H_[0,6]
                H_[0,7] = 2.0*DY1*DY8
                H_[7,0] = H_[0,7]
                H_[0,8] = 2.0*DY1*DY9
                H_[8,0] = H_[0,8]
                H_[0,9] = 2.0*DY1*DY10
                H_[9,0] = H_[0,9]
                H_[0,10] = 2.0*DY1*DY11
                H_[10,0] = H_[0,10]
                H_[0,11] = 2.0*DY1*DY12
                H_[11,0] = H_[0,11]
                H_[0,12] = 2.0*DY1*DY13
                H_[12,0] = H_[0,12]
                H_[0,13] = 2.0*DY1*DY14
                H_[13,0] = H_[0,13]
                H_[0,14] = 2.0*DY1*DY15
                H_[14,0] = H_[0,14]
                H_[0,15] = 2.0*DY1*DY16
                H_[15,0] = H_[0,15]
                H_[0,16] = 2.0*DY1*DY17
                H_[16,0] = H_[0,16]
                H_[0,17] = 2.0*DY1*DY18
                H_[17,0] = H_[0,17]
                H_[0,18] = 2.0*DY1*DY19
                H_[18,0] = H_[0,18]
                H_[0,19] = 2.0*DY1*DY20
                H_[19,0] = H_[0,19]
                H_[1,1] = 2.0*DY2*DY2+TWOSUM*D2Y2
                H_[1,2] = 2.0*DY2*DY3
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0*DY2*DY4
                H_[3,1] = H_[1,3]
                H_[1,4] = 2.0*DY2*DY5
                H_[4,1] = H_[1,4]
                H_[1,5] = 2.0*DY2*DY6
                H_[5,1] = H_[1,5]
                H_[1,6] = 2.0*DY2*DY7
                H_[6,1] = H_[1,6]
                H_[1,7] = 2.0*DY2*DY8
                H_[7,1] = H_[1,7]
                H_[1,8] = 2.0*DY2*DY9
                H_[8,1] = H_[1,8]
                H_[1,9] = 2.0*DY2*DY10
                H_[9,1] = H_[1,9]
                H_[1,10] = 2.0*DY2*DY11
                H_[10,1] = H_[1,10]
                H_[1,11] = 2.0*DY2*DY12
                H_[11,1] = H_[1,11]
                H_[1,12] = 2.0*DY2*DY13
                H_[12,1] = H_[1,12]
                H_[1,13] = 2.0*DY2*DY14
                H_[13,1] = H_[1,13]
                H_[1,14] = 2.0*DY2*DY15
                H_[14,1] = H_[1,14]
                H_[1,15] = 2.0*DY2*DY16
                H_[15,1] = H_[1,15]
                H_[1,16] = 2.0*DY2*DY17
                H_[16,1] = H_[1,16]
                H_[1,17] = 2.0*DY2*DY18
                H_[17,1] = H_[1,17]
                H_[1,18] = 2.0*DY2*DY19
                H_[18,1] = H_[1,18]
                H_[1,19] = 2.0*DY2*DY20
                H_[19,1] = H_[1,19]
                H_[2,2] = 2.0*DY3*DY3+TWOSUM*D2Y3
                H_[2,3] = 2.0*DY3*DY4
                H_[3,2] = H_[2,3]
                H_[2,4] = 2.0*DY3*DY5
                H_[4,2] = H_[2,4]
                H_[2,5] = 2.0*DY3*DY6
                H_[5,2] = H_[2,5]
                H_[2,6] = 2.0*DY3*DY7
                H_[6,2] = H_[2,6]
                H_[2,7] = 2.0*DY3*DY8
                H_[7,2] = H_[2,7]
                H_[2,8] = 2.0*DY3*DY9
                H_[8,2] = H_[2,8]
                H_[2,9] = 2.0*DY3*DY10
                H_[9,2] = H_[2,9]
                H_[2,10] = 2.0*DY3*DY11
                H_[10,2] = H_[2,10]
                H_[2,11] = 2.0*DY3*DY12
                H_[11,2] = H_[2,11]
                H_[2,12] = 2.0*DY3*DY13
                H_[12,2] = H_[2,12]
                H_[2,13] = 2.0*DY3*DY14
                H_[13,2] = H_[2,13]
                H_[2,14] = 2.0*DY3*DY15
                H_[14,2] = H_[2,14]
                H_[2,15] = 2.0*DY3*DY16
                H_[15,2] = H_[2,15]
                H_[2,16] = 2.0*DY3*DY17
                H_[16,2] = H_[2,16]
                H_[2,17] = 2.0*DY3*DY18
                H_[17,2] = H_[2,17]
                H_[2,18] = 2.0*DY3*DY19
                H_[18,2] = H_[2,18]
                H_[2,19] = 2.0*DY3*DY20
                H_[19,2] = H_[2,19]
                H_[3,3] = 2.0*DY4*DY4+TWOSUM*D2Y4
                H_[3,4] = 2.0*DY4*DY5
                H_[4,3] = H_[3,4]
                H_[3,5] = 2.0*DY4*DY6
                H_[5,3] = H_[3,5]
                H_[3,6] = 2.0*DY4*DY7
                H_[6,3] = H_[3,6]
                H_[3,7] = 2.0*DY4*DY8
                H_[7,3] = H_[3,7]
                H_[3,8] = 2.0*DY4*DY9
                H_[8,3] = H_[3,8]
                H_[3,9] = 2.0*DY4*DY10
                H_[9,3] = H_[3,9]
                H_[3,10] = 2.0*DY4*DY11
                H_[10,3] = H_[3,10]
                H_[3,11] = 2.0*DY4*DY12
                H_[11,3] = H_[3,11]
                H_[3,12] = 2.0*DY4*DY13
                H_[12,3] = H_[3,12]
                H_[3,13] = 2.0*DY4*DY14
                H_[13,3] = H_[3,13]
                H_[3,14] = 2.0*DY4*DY15
                H_[14,3] = H_[3,14]
                H_[3,15] = 2.0*DY4*DY16
                H_[15,3] = H_[3,15]
                H_[3,16] = 2.0*DY4*DY17
                H_[16,3] = H_[3,16]
                H_[3,17] = 2.0*DY4*DY18
                H_[17,3] = H_[3,17]
                H_[3,18] = 2.0*DY4*DY19
                H_[18,3] = H_[3,18]
                H_[3,19] = 2.0*DY4*DY20
                H_[19,3] = H_[3,19]
                H_[4,4] = 2.0*DY5*DY5+TWOSUM*D2Y5
                H_[4,5] = 2.0*DY5*DY6
                H_[5,4] = H_[4,5]
                H_[4,6] = 2.0*DY5*DY7
                H_[6,4] = H_[4,6]
                H_[4,7] = 2.0*DY5*DY8
                H_[7,4] = H_[4,7]
                H_[4,8] = 2.0*DY5*DY9
                H_[8,4] = H_[4,8]
                H_[4,9] = 2.0*DY5*DY10
                H_[9,4] = H_[4,9]
                H_[4,10] = 2.0*DY5*DY11
                H_[10,4] = H_[4,10]
                H_[4,11] = 2.0*DY5*DY12
                H_[11,4] = H_[4,11]
                H_[4,12] = 2.0*DY5*DY13
                H_[12,4] = H_[4,12]
                H_[4,13] = 2.0*DY5*DY14
                H_[13,4] = H_[4,13]
                H_[4,14] = 2.0*DY5*DY15
                H_[14,4] = H_[4,14]
                H_[4,15] = 2.0*DY5*DY16
                H_[15,4] = H_[4,15]
                H_[4,16] = 2.0*DY5*DY17
                H_[16,4] = H_[4,16]
                H_[4,17] = 2.0*DY5*DY18
                H_[17,4] = H_[4,17]
                H_[4,18] = 2.0*DY5*DY19
                H_[18,4] = H_[4,18]
                H_[4,19] = 2.0*DY5*DY20
                H_[19,4] = H_[4,19]
                H_[5,5] = 2.0*DY6*DY6+TWOSUM*D2Y6
                H_[5,6] = 2.0*DY6*DY7
                H_[6,5] = H_[5,6]
                H_[5,7] = 2.0*DY6*DY8
                H_[7,5] = H_[5,7]
                H_[5,8] = 2.0*DY6*DY9
                H_[8,5] = H_[5,8]
                H_[5,9] = 2.0*DY6*DY10
                H_[9,5] = H_[5,9]
                H_[5,10] = 2.0*DY6*DY11
                H_[10,5] = H_[5,10]
                H_[5,11] = 2.0*DY6*DY12
                H_[11,5] = H_[5,11]
                H_[5,12] = 2.0*DY6*DY13
                H_[12,5] = H_[5,12]
                H_[5,13] = 2.0*DY6*DY14
                H_[13,5] = H_[5,13]
                H_[5,14] = 2.0*DY6*DY15
                H_[14,5] = H_[5,14]
                H_[5,15] = 2.0*DY6*DY16
                H_[15,5] = H_[5,15]
                H_[5,16] = 2.0*DY6*DY17
                H_[16,5] = H_[5,16]
                H_[5,17] = 2.0*DY6*DY18
                H_[17,5] = H_[5,17]
                H_[5,18] = 2.0*DY6*DY19
                H_[18,5] = H_[5,18]
                H_[5,19] = 2.0*DY6*DY20
                H_[19,5] = H_[5,19]
                H_[6,6] = 2.0*DY7*DY7+TWOSUM*D2Y7
                H_[6,7] = 2.0*DY7*DY8
                H_[7,6] = H_[6,7]
                H_[6,8] = 2.0*DY7*DY9
                H_[8,6] = H_[6,8]
                H_[6,9] = 2.0*DY7*DY10
                H_[9,6] = H_[6,9]
                H_[6,10] = 2.0*DY7*DY11
                H_[10,6] = H_[6,10]
                H_[6,11] = 2.0*DY7*DY12
                H_[11,6] = H_[6,11]
                H_[6,12] = 2.0*DY7*DY13
                H_[12,6] = H_[6,12]
                H_[6,13] = 2.0*DY7*DY14
                H_[13,6] = H_[6,13]
                H_[6,14] = 2.0*DY7*DY15
                H_[14,6] = H_[6,14]
                H_[6,15] = 2.0*DY7*DY16
                H_[15,6] = H_[6,15]
                H_[6,16] = 2.0*DY7*DY17
                H_[16,6] = H_[6,16]
                H_[6,17] = 2.0*DY7*DY18
                H_[17,6] = H_[6,17]
                H_[6,18] = 2.0*DY7*DY19
                H_[18,6] = H_[6,18]
                H_[6,19] = 2.0*DY7*DY20
                H_[19,6] = H_[6,19]
                H_[7,7] = 2.0*DY8*DY8+TWOSUM*D2Y8
                H_[7,8] = 2.0*DY8*DY9
                H_[8,7] = H_[7,8]
                H_[7,9] = 2.0*DY8*DY10
                H_[9,7] = H_[7,9]
                H_[7,10] = 2.0*DY8*DY11
                H_[10,7] = H_[7,10]
                H_[7,11] = 2.0*DY8*DY12
                H_[11,7] = H_[7,11]
                H_[7,12] = 2.0*DY8*DY13
                H_[12,7] = H_[7,12]
                H_[7,13] = 2.0*DY8*DY14
                H_[13,7] = H_[7,13]
                H_[7,14] = 2.0*DY8*DY15
                H_[14,7] = H_[7,14]
                H_[7,15] = 2.0*DY8*DY16
                H_[15,7] = H_[7,15]
                H_[7,16] = 2.0*DY8*DY17
                H_[16,7] = H_[7,16]
                H_[7,17] = 2.0*DY8*DY18
                H_[17,7] = H_[7,17]
                H_[7,18] = 2.0*DY8*DY19
                H_[18,7] = H_[7,18]
                H_[7,19] = 2.0*DY8*DY20
                H_[19,7] = H_[7,19]
                H_[8,8] = 2.0*DY9*DY9+TWOSUM*D2Y9
                H_[8,9] = 2.0*DY9*DY10
                H_[9,8] = H_[8,9]
                H_[8,10] = 2.0*DY9*DY11
                H_[10,8] = H_[8,10]
                H_[8,11] = 2.0*DY9*DY12
                H_[11,8] = H_[8,11]
                H_[8,12] = 2.0*DY9*DY13
                H_[12,8] = H_[8,12]
                H_[8,13] = 2.0*DY9*DY14
                H_[13,8] = H_[8,13]
                H_[8,14] = 2.0*DY9*DY15
                H_[14,8] = H_[8,14]
                H_[8,15] = 2.0*DY9*DY16
                H_[15,8] = H_[8,15]
                H_[8,16] = 2.0*DY9*DY17
                H_[16,8] = H_[8,16]
                H_[8,17] = 2.0*DY9*DY18
                H_[17,8] = H_[8,17]
                H_[8,18] = 2.0*DY9*DY19
                H_[18,8] = H_[8,18]
                H_[8,19] = 2.0*DY9*DY20
                H_[19,8] = H_[8,19]
                H_[9,9] = 2.0*DY10*DY10+TWOSUM*D2Y10
                H_[9,10] = 2.0*DY10*DY11
                H_[10,9] = H_[9,10]
                H_[9,11] = 2.0*DY10*DY12
                H_[11,9] = H_[9,11]
                H_[9,12] = 2.0*DY10*DY13
                H_[12,9] = H_[9,12]
                H_[9,13] = 2.0*DY10*DY14
                H_[13,9] = H_[9,13]
                H_[9,14] = 2.0*DY10*DY15
                H_[14,9] = H_[9,14]
                H_[9,15] = 2.0*DY10*DY16
                H_[15,9] = H_[9,15]
                H_[9,16] = 2.0*DY10*DY17
                H_[16,9] = H_[9,16]
                H_[9,17] = 2.0*DY10*DY18
                H_[17,9] = H_[9,17]
                H_[9,18] = 2.0*DY10*DY19
                H_[18,9] = H_[9,18]
                H_[9,19] = 2.0*DY10*DY20
                H_[19,9] = H_[9,19]
                H_[10,10] = 2.0*DY11*DY11+TWOSUM*D2Y11
                H_[10,11] = 2.0*DY11*DY12
                H_[11,10] = H_[10,11]
                H_[10,12] = 2.0*DY11*DY13
                H_[12,10] = H_[10,12]
                H_[10,13] = 2.0*DY11*DY14
                H_[13,10] = H_[10,13]
                H_[10,14] = 2.0*DY11*DY15
                H_[14,10] = H_[10,14]
                H_[10,15] = 2.0*DY11*DY16
                H_[15,10] = H_[10,15]
                H_[10,16] = 2.0*DY11*DY17
                H_[16,10] = H_[10,16]
                H_[10,17] = 2.0*DY11*DY18
                H_[17,10] = H_[10,17]
                H_[10,18] = 2.0*DY11*DY19
                H_[18,10] = H_[10,18]
                H_[10,19] = 2.0*DY11*DY20
                H_[19,10] = H_[10,19]
                H_[11,11] = 2.0*DY12*DY12+TWOSUM*D2Y12
                H_[11,12] = 2.0*DY12*DY13
                H_[12,11] = H_[11,12]
                H_[11,13] = 2.0*DY12*DY14
                H_[13,11] = H_[11,13]
                H_[11,14] = 2.0*DY12*DY15
                H_[14,11] = H_[11,14]
                H_[11,15] = 2.0*DY12*DY16
                H_[15,11] = H_[11,15]
                H_[11,16] = 2.0*DY12*DY17
                H_[16,11] = H_[11,16]
                H_[11,17] = 2.0*DY12*DY18
                H_[17,11] = H_[11,17]
                H_[11,18] = 2.0*DY12*DY19
                H_[18,11] = H_[11,18]
                H_[11,19] = 2.0*DY12*DY20
                H_[19,11] = H_[11,19]
                H_[12,12] = 2.0*DY13*DY13+TWOSUM*D2Y13
                H_[12,13] = 2.0*DY13*DY14
                H_[13,12] = H_[12,13]
                H_[12,14] = 2.0*DY13*DY15
                H_[14,12] = H_[12,14]
                H_[12,15] = 2.0*DY13*DY16
                H_[15,12] = H_[12,15]
                H_[12,16] = 2.0*DY13*DY17
                H_[16,12] = H_[12,16]
                H_[12,17] = 2.0*DY13*DY18
                H_[17,12] = H_[12,17]
                H_[12,18] = 2.0*DY13*DY19
                H_[18,12] = H_[12,18]
                H_[12,19] = 2.0*DY13*DY20
                H_[19,12] = H_[12,19]
                H_[13,13] = 2.0*DY14*DY14+TWOSUM*D2Y14
                H_[13,14] = 2.0*DY14*DY15
                H_[14,13] = H_[13,14]
                H_[13,15] = 2.0*DY14*DY16
                H_[15,13] = H_[13,15]
                H_[13,16] = 2.0*DY14*DY17
                H_[16,13] = H_[13,16]
                H_[13,17] = 2.0*DY14*DY18
                H_[17,13] = H_[13,17]
                H_[13,18] = 2.0*DY14*DY19
                H_[18,13] = H_[13,18]
                H_[13,19] = 2.0*DY14*DY20
                H_[19,13] = H_[13,19]
                H_[14,14] = 2.0*DY15*DY15+TWOSUM*D2Y15
                H_[14,15] = 2.0*DY15*DY16
                H_[15,14] = H_[14,15]
                H_[14,16] = 2.0*DY15*DY17
                H_[16,14] = H_[14,16]
                H_[14,17] = 2.0*DY15*DY18
                H_[17,14] = H_[14,17]
                H_[14,18] = 2.0*DY15*DY19
                H_[18,14] = H_[14,18]
                H_[14,19] = 2.0*DY15*DY20
                H_[19,14] = H_[14,19]
                H_[15,15] = 2.0*DY16*DY16+TWOSUM*D2Y16
                H_[15,16] = 2.0*DY16*DY17
                H_[16,15] = H_[15,16]
                H_[15,17] = 2.0*DY16*DY18
                H_[17,15] = H_[15,17]
                H_[15,18] = 2.0*DY16*DY19
                H_[18,15] = H_[15,18]
                H_[15,19] = 2.0*DY16*DY20
                H_[19,15] = H_[15,19]
                H_[16,16] = 2.0*DY17*DY17+TWOSUM*D2Y17
                H_[16,17] = 2.0*DY17*DY18
                H_[17,16] = H_[16,17]
                H_[16,18] = 2.0*DY17*DY19
                H_[18,16] = H_[16,18]
                H_[16,19] = 2.0*DY17*DY20
                H_[19,16] = H_[16,19]
                H_[17,17] = 2.0*DY18*DY18+TWOSUM*D2Y18
                H_[17,18] = 2.0*DY18*DY19
                H_[18,17] = H_[17,18]
                H_[17,19] = 2.0*DY18*DY20
                H_[19,17] = H_[17,19]
                H_[18,18] = 2.0*DY19*DY19+TWOSUM*D2Y19
                H_[18,19] = 2.0*DY19*DY20
                H_[19,18] = H_[18,19]
                H_[19,19] = 2.0*DY20*DY20+TWOSUM*D2Y20
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eQR(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]**4
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 4.0*EV_[0]**3
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 12.0*EV_[0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def en3P(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[1]*EV_[2]+2.0*EV_[2]*EV_[2]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[1]*EV_[2]
            g_[1] = EV_[0]*EV_[2]
            g_[2] = EV_[0]*EV_[1]+4.0*EV_[2]
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = EV_[2]
                H_[1,0] = H_[0,1]
                H_[0,2] = EV_[1]
                H_[2,0] = H_[0,2]
                H_[1,2] = EV_[0]
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

