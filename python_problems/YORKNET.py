from s2mpjlib import *
class  YORKNET(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A problem arising in the modelling of the Yorkshire water system.
# 
#    Source:
#    an problem submitted for the LANCELOT licence.
# 
#    SIF input: B. Ulanicki, Water Software Systems,De Montfort University,
#               The Gateway, Leicester LE1 9BH, UK.
#               e-mail: bul@uk.ac.dmu * Tel no.0533 577070
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-CSOR2-AY-312-256"
# 
# DECLARE CONSTANTS DESCRIBING NETWORK
# 
# STANDARD DECLARATIONS
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'YORKNET'

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
        v_['2'] = 2
        v_['3'] = 3
        v_['4'] = 4
        v_['5'] = 5
        v_['6'] = 6
        v_['7'] = 7
        v_['8'] = 8
        v_['9'] = 9
        v_['10'] = 10
        v_['11'] = 11
        v_['12'] = 12
        v_['13'] = 13
        v_['14'] = 14
        v_['15'] = 15
        v_['16'] = 16
        v_['17'] = 17
        v_['18'] = 18
        v_['19'] = 19
        v_['ONE'] = 1
        v_['NSTEP'] = 8
        v_['NSTEP+1'] = 8+v_['NSTEP']
        v_['ST1'] = 0.125
        v_['ST2'] = 0.125
        v_['ST3'] = 0.125
        v_['ST4'] = 0.125
        v_['ST5'] = 0.125
        v_['ST6'] = 0.125
        v_['ST7'] = 0.125
        v_['ST8'] = 0.125
        v_['MST1'] = 0.0-v_['ST1']
        v_['MST2'] = 0.0-v_['ST2']
        v_['MST3'] = 0.0-v_['ST3']
        v_['MST4'] = 0.0-v_['ST4']
        v_['MST5'] = 0.0-v_['ST5']
        v_['MST6'] = 0.0-v_['ST6']
        v_['MST7'] = 0.0-v_['ST7']
        v_['MST8'] = 0.0-v_['ST8']
        v_['TFL1'] = 1.87
        v_['TFH1'] = 4.03
        v_['TF1'] = 0.0+v_['TFL1']
        v_['TF2'] = 0.0+v_['TFL1']
        v_['TF3'] = 0.0+v_['TFH1']
        v_['TF4'] = 0.0+v_['TFH1']
        v_['TF5'] = 0.0+v_['TFH1']
        v_['TF6'] = 0.0+v_['TFH1']
        v_['TF7'] = 0.0+v_['TFH1']
        v_['TF8'] = 0.0+v_['TFH1']
        v_['BP1'] = 24.0*v_['ST1']
        v_['BP2'] = 24.0*v_['ST2']
        v_['BP3'] = 24.0*v_['ST3']
        v_['BP4'] = 24.0*v_['ST4']
        v_['BP5'] = 24.0*v_['ST5']
        v_['BP6'] = 24.0*v_['ST6']
        v_['BP7'] = 24.0*v_['ST7']
        v_['BP8'] = 24.0*v_['ST8']
        v_['CP1'] = v_['BP1']*v_['TF1']
        v_['CP2'] = v_['BP2']*v_['TF2']
        v_['CP3'] = v_['BP3']*v_['TF3']
        v_['CP4'] = v_['BP4']*v_['TF4']
        v_['CP5'] = v_['BP5']*v_['TF5']
        v_['CP6'] = v_['BP6']*v_['TF6']
        v_['CP7'] = v_['BP7']*v_['TF7']
        v_['CP8'] = v_['BP8']*v_['TF8']
        v_['NEL'] = 15
        v_['SuVAL'] = 9
        v_['EuVAL'] = 10
        v_['SuPIPE'] = 1
        v_['EuPIPE'] = 8
        v_['MGINV1'] = -3.365170e-3
        v_['MGINV2'] = -2.314284e-2
        v_['MGINV3'] = -6.631203e-3
        v_['MGINV4'] = -1.702093e-3
        v_['MGINV5'] = -1.205983e-2
        v_['MGINV6'] = -9.617776e-4
        v_['MGINV7'] = -1.392046e-5
        v_['MGINV8'] = -4.411625e-3
        v_['MGINV9'] = -2.019250e-3
        v_['MGINV10'] = -2.288437e-3
        v_['NPMP'] = 5
        v_['SuPMP'] = 11
        v_['EuPMP'] = 15
        v_['CONA11'] = -.035520
        v_['CONB11'] = -.054720
        v_['CONC11'] = 99.80
        v_['CONA12'] = -0.07475
        v_['CONB12'] = -9.05
        v_['CONC12'] = 110
        v_['CONA13'] = -.042420
        v_['CONB13'] = -.005370
        v_['CONC13'] = 175.29
        v_['CONA14'] = -.040733
        v_['CONB14'] = -.032036
        v_['CONC14'] = 139.6
        v_['CONA15'] = -.167495
        v_['CONB15'] = -.0019
        v_['CONC15'] = 139.6
        v_['NND'] = 13
        v_['D1'] = 0.0
        v_['D4'] = -33.0
        v_['D2'] = 0.0
        v_['D3'] = -55.0
        v_['D5'] = 0.0
        v_['D6'] = 0.0
        v_['D7'] = 0.0
        v_['D8'] = -25.0
        v_['D9'] = 0.0
        v_['D10'] = -17.0
        v_['D11'] = 0.0
        v_['D12'] = 0.0
        v_['D13'] = 0.0
        v_['SuRES'] = 1
        v_['EuRES'] = 4
        v_['EuRES+1'] = 1+v_['EuRES']
        v_['HGT1'] = 5.77
        v_['HGT2'] = 3.00
        v_['HGT3'] = 131.08
        v_['HGT4'] = 44.0
        v_['MXHGT1'] = 9.60
        v_['MXHGT2'] = 7.89
        v_['MXHGT3'] = 138.76
        v_['MXHGT4'] = 53.34
        v_['XAR1'] = 1.599
        v_['XAR2'] = 4.6421
        v_['XAR3'] = 30.2307
        v_['XAR4'] = 5.3938
        v_['MXAR1'] = 0.0-v_['XAR1']
        v_['MXAR2'] = 0.0-v_['XAR2']
        v_['MXAR3'] = 0.0-v_['XAR3']
        v_['MXAR4'] = 0.0-v_['XAR4']
        v_['RXAR1'] = 1.0/v_['XAR1']
        v_['MRXAR1'] = 0.0-v_['RXAR1']
        v_['RXAR2'] = 1.0/v_['XAR2']
        v_['MRXAR2'] = 0.0-v_['RXAR2']
        v_['RXAR3'] = 1.0/v_['XAR3']
        v_['MRXAR3'] = 0.0-v_['RXAR3']
        v_['RXAR4'] = 1.0/v_['XAR4']
        v_['MRXAR4'] = 0.0-v_['RXAR4']
        v_['HTXAR1'] = v_['XAR1']*v_['HGT1']
        v_['HTXAR2'] = v_['XAR2']*v_['HGT2']
        v_['HTXAR3'] = v_['XAR3']*v_['HGT3']
        v_['HTXAR4'] = v_['XAR4']*v_['HGT4']
        v_['STHGT1'] = 8.5
        v_['STHGT2'] = 6.0
        v_['STHGT3'] = 135.6
        v_['STHGT4'] = 48.5
        v_['STVOL1'] = v_['STHGT1']*v_['XAR1']
        v_['STVOL2'] = v_['STHGT2']*v_['XAR2']
        v_['STVOL3'] = v_['STHGT3']*v_['XAR3']
        v_['STVOL4'] = v_['STHGT4']*v_['XAR4']
        v_['MSTVOL1'] = 0.0-v_['STVOL1']
        v_['MSTVOL2'] = 0.0-v_['STVOL2']
        v_['MSTVOL3'] = 0.0-v_['STVOL3']
        v_['MSTVOL4'] = 0.0-v_['STVOL4']
        v_['WMN1'] = 3.764
        v_['WMN2'] = 11.35
        v_['WMN3'] = 156.648
        v_['WMN4'] = 45.929
        v_['WMX1'] = 5.646
        v_['WMX2'] = 22.133
        v_['WMX3'] = 223.489
        v_['WMX4'] = 61.876
        v_['H1'] = 8.99
        v_['H2'] = 52.84
        v_['H3'] = 138.31
        v_['H4'] = 5.67
        v_['W1'] = 4.728
        v_['W2'] = 15.601
        v_['W3'] = 190.648
        v_['W4'] = 55.00
        v_['SuTW'] = 1
        v_['EuTW'] = 2
        v_['BCTW1'] = 28.34
        v_['BCTW2'] = 18.86
        v_['BCTW3'] = 1.0
        v_['BCTW4'] = 1.0
        v_['BCTW5'] = 28.34
        v_['BCTW6'] = 18.86
        v_['BCTW7'] = 1.0
        v_['BCTW8'] = 1.0
        v_['CTW1,1'] = v_['BCTW1']*v_['ST1']
        v_['CTW1,2'] = v_['BCTW1']*v_['ST2']
        v_['CTW1,3'] = v_['BCTW1']*v_['ST3']
        v_['CTW1,4'] = v_['BCTW1']*v_['ST4']
        v_['CTW1,5'] = v_['BCTW1']*v_['ST5']
        v_['CTW1,6'] = v_['BCTW1']*v_['ST6']
        v_['CTW1,7'] = v_['BCTW1']*v_['ST7']
        v_['CTW1,8'] = v_['BCTW1']*v_['ST8']
        v_['CTW2,1'] = v_['BCTW2']*v_['ST1']
        v_['CTW2,2'] = v_['BCTW2']*v_['ST2']
        v_['CTW2,3'] = v_['BCTW2']*v_['ST3']
        v_['CTW2,4'] = v_['BCTW2']*v_['ST4']
        v_['CTW2,5'] = v_['BCTW2']*v_['ST5']
        v_['CTW2,6'] = v_['BCTW2']*v_['ST6']
        v_['CTW2,7'] = v_['BCTW2']*v_['ST7']
        v_['CTW2,8'] = v_['BCTW2']*v_['ST8']
        v_['WS1'] = 94.67
        v_['WS2'] = 30.87
        v_['QSMX1'] = 500.0
        v_['QSMX2'] = 40.0
        v_['SCuQ'] = 1.0
        v_['SCuH'] = 1.0
        v_['SCuV'] = 1.0
        v_['PIPEuSC'] = 1.0
        v_['VALuSC'] = 1.0
        v_['NDuSC'] = 1.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        irA          = np.array([],dtype=int)
        icA          = np.array([],dtype=int)
        valA         = np.array([],dtype=float)
        for I in range(int(v_['ONE']),int(v_['NSTEP'])+1):
            for J in range(int(v_['ONE']),int(v_['NEL'])+1):
                [iv,ix_,_] = s2mpj_ii('q'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'q'+str(J)+','+str(I))
            for J in range(int(v_['ONE']),int(v_['NND'])+1):
                [iv,ix_,_] = s2mpj_ii('h'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'h'+str(J)+','+str(I))
            for J in range(int(v_['SuRES']),int(v_['EuRES'])+1):
                [iv,ix_,_] = s2mpj_ii('qr'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'qr'+str(J)+','+str(I))
            for J in range(int(v_['SuTW']),int(v_['EuTW'])+1):
                [iv,ix_,_] = s2mpj_ii('qs'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'qs'+str(J)+','+str(I))
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                [iv,ix_,_] = s2mpj_ii('u'+str(J)+','+str(I),ix_)
                self.xnames=arrset(self.xnames,iv,'u'+str(J)+','+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames       = np.array([])
        self.cnames  = np.array([])
        gtype        = np.array([])
        for I in range(int(v_['ONE']),int(v_['NSTEP'])+1):
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qr'+str(int(v_['1']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qs'+str(int(v_['1']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['13']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['2']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['2']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qr'+str(int(v_['2']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qs'+str(int(v_['2']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['2']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['2']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['14']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['15']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['3']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['3']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qr'+str(int(v_['3']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['5']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['3']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['3']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['6']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['4']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['4']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qr'+str(int(v_['4']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['12']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['4']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['4']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['8']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['9']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['5']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['5']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['13']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['1']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['6']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['6']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['2']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['3']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['6']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['6']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['12']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['7']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['7']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['2']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['7']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['7']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['7']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['11']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['8']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['8']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['3']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['7']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['9']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['9']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['4']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['5']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['9']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['9']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['11']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['10']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['10']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['4']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['6']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['11']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['11']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['10']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['14']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['11']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['11']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['15']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['12']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['12']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['8']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['10']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            [ig,ig_,_] = s2mpj_ii('ND'+str(int(v_['13']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'ND'+str(int(v_['13']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['9']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['q'+str(int(v_['1']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
        for I in range(int(v_['1']),int(v_['NSTEP'])+1):
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['1']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['1']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['13']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['5']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['2']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['2']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['7']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['6']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['3']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['3']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['8']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['6']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['4']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['4']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['10']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['9']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['5']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['5']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['3']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['9']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['6']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['6']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['10']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['3']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['7']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['7']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['7']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['8']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['8']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['8']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['4']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['12']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['9']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['9']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['4']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['13']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            [ig,ig_,_] = s2mpj_ii('EL'+str(int(v_['10']))+','+str(I),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'EL'+str(int(v_['10']))+','+str(I))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['12']))+','+str(I)]])
            valA = np.append(valA,float(-1.0))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(int(v_['11']))+','+str(I)]])
            valA = np.append(valA,float(1.0))
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                [ig,ig_,_] = s2mpj_ii('EL'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'EL'+str(J)+','+str(I))
        for I in range(int(v_['ONE']),int(v_['NSTEP'])+1):
            for J in range(int(v_['SuTW']),int(v_['EuTW'])+1):
                [ig,ig_,_] = s2mpj_ii('TWC'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['qs'+str(J)+','+str(I)]])
                valA = np.append(valA,float(v_['CTW'+str(J)+','+str(I)]))
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                [ig,ig_,_] = s2mpj_ii('PC'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'<>')
        for J in range(int(v_['SuRES']),int(v_['EuRES'])+1):
            [ig,ig_,_] = s2mpj_ii('RD'+str(J)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'RD'+str(J)+','+str(int(v_['1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['h'+str(J)+','+str(int(v_['1']))]])
            valA = np.append(valA,float(v_['MXAR'+str(J)]))
            [ig,ig_,_] = s2mpj_ii('RD'+str(J)+','+str(int(v_['1'])),ig_)
            gtype = arrset(gtype,ig,'==')
            cnames = arrset(cnames,ig,'RD'+str(J)+','+str(int(v_['1'])))
            irA  = np.append(irA,[ig])
            icA  = np.append(icA,[ix_['qr'+str(J)+','+str(int(v_['1']))]])
            valA = np.append(valA,float(v_['MST'+str(int(v_['1']))]))
        for I in range(int(v_['2']),int(v_['NSTEP'])+1):
            for J in range(int(v_['SuRES']),int(v_['EuRES'])+1):
                v_['A'] = -1+I
                [ig,ig_,_] = s2mpj_ii('RD'+str(J)+','+str(I),ig_)
                gtype = arrset(gtype,ig,'==')
                cnames = arrset(cnames,ig,'RD'+str(J)+','+str(I))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['h'+str(J)+','+str(I)]])
                valA = np.append(valA,float(v_['MXAR'+str(J)]))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['h'+str(J)+','+str(int(v_['A']))]])
                valA = np.append(valA,float(v_['XAR'+str(J)]))
                irA  = np.append(irA,[ig])
                icA  = np.append(icA,[ix_['qr'+str(J)+','+str(I)]])
                valA = np.append(valA,float(v_['MST'+str(I)]))
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
        for I in range(int(v_['ONE']),int(v_['NSTEP'])+1):
            for J in range(int(v_['ONE']),int(v_['NND'])+1):
                self.gconst  = (
                      arrset(self.gconst,ig_['ND'+str(J)+','+str(I)],float(v_['D'+str(J)])))
        for J in range(int(v_['SuRES']),int(v_['EuRES'])+1):
            self.gconst  = (
                  arrset(self.gconst,ig_['RD'+str(J)+','+str(int(v_['1']))],float(v_['MSTVOL'+str(J)])))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.zeros((self.n,1))
        self.xupper = np.full((self.n,1),float('inf'))
        for I in range(int(v_['ONE']),int(v_['NSTEP'])+1):
            for J in range(int(v_['SuRES']),int(v_['EuRES'])+1):
                self.xlower[ix_['h'+str(J)+','+str(I)]] = v_['HGT'+str(J)]
                self.xupper[ix_['h'+str(J)+','+str(I)]] = v_['MXHGT'+str(J)]
                self.xlower[ix_['qr'+str(J)+','+str(I)]] = -float('Inf')
                self.xupper[ix_['qr'+str(J)+','+str(I)]] = +float('Inf')
            for J in range(int(v_['EuRES+1']),int(v_['NND'])+1):
                self.xlower[ix_['h'+str(J)+','+str(I)]] = 0.0
            for J in range(int(v_['SuTW']),int(v_['EuTW'])+1):
                self.xlower[ix_['qs'+str(J)+','+str(I)]] = 0.0
                self.xupper[ix_['qs'+str(J)+','+str(I)]] = v_['QSMX'+str(J)]
            for J in range(int(v_['ONE']),int(v_['NEL'])+1):
                self.xlower[ix_['q'+str(J)+','+str(I)]] = -float('Inf')
                self.xupper[ix_['q'+str(J)+','+str(I)]] = +float('Inf')
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                self.xlower[ix_['u'+str(J)+','+str(I)]] = 0.0
                self.xupper[ix_['u'+str(J)+','+str(I)]] = 7.0
                self.xlower[ix_['q'+str(J)+','+str(I)]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.full((self.n,1),float(0.0))
        self.y0 = np.full((self.m,1),float(0.0))
        for I in range(int(v_['1']),int(v_['NSTEP'])+1):
            for J in range(int(v_['EuRES+1']),int(v_['NND'])+1):
                self.x0[ix_['h'+str(J)+','+str(I)]] = float(0.0)
            for J in range(int(v_['ONE']),int(v_['NEL'])+1):
                self.x0[ix_['q'+str(J)+','+str(I)]] = float(20.0)
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                self.x0[ix_['u'+str(J)+','+str(I)]] = float(3.0)
            self.x0[ix_['qs'+str(int(v_['1']))+','+str(I)]] = float(25.0)
            self.x0[ix_['qs'+str(int(v_['2']))+','+str(I)]] = float(50.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eXSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eXZ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Z')
        [it,iet_,_] = s2mpj_ii( 'eXPOW', iet_)
        elftv = loaset(elftv,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eXMYZSQ', iet_)
        elftv = loaset(elftv,it,0,'X')
        elftv = loaset(elftv,it,1,'Y')
        elftv = loaset(elftv,it,2,'Z')
        [it,iet_,_] = s2mpj_ii( 'ePMP1', iet_)
        elftv = loaset(elftv,it,0,'Q')
        [it,iet_,_] = s2mpj_ii( 'ePMP2', iet_)
        elftv = loaset(elftv,it,0,'Q')
        [it,iet_,_] = s2mpj_ii( 'ePMP3', iet_)
        elftv = loaset(elftv,it,0,'Q')
        [it,iet_,_] = s2mpj_ii( 'ePMP4', iet_)
        elftv = loaset(elftv,it,0,'Q')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        for I in range(int(v_['1']),int(v_['NSTEP'])+1):
            for J in range(int(v_['SuPIPE']),int(v_['EuPIPE'])+1):
                ename = 'EA'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXPOW')
                ielftype = arrset(ielftype,ie,iet_["eXPOW"])
                vname = 'q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['SuVAL']),int(v_['EuVAL'])+1):
                ename = 'EA'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXPOW')
                ielftype = arrset(ielftype,ie,iet_["eXPOW"])
                vname = 'q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                ename = 'EA'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXSQ')
                ielftype = arrset(ielftype,ie,iet_["eXSQ"])
                vname = 'q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'EB'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXZ')
                ielftype = arrset(ielftype,ie,iet_["eXZ"])
                vname = 'q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                vname = 'u'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='Z')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'EC'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXSQ')
                ielftype = arrset(ielftype,ie,iet_["eXSQ"])
                vname = 'u'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='X')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
                ename = 'EH'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                self.elftype = arrset(self.elftype,ie,'eXMYZSQ')
                ielftype = arrset(ielftype,ie,iet_["eXMYZSQ"])
                vname = 'u'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='Z')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['11']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['7']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['11']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['9']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['12']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['4']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['12']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['6']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['13']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['1']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['13']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['5']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['14']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['2']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['14']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['11']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['15']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['2']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='X')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            ename = 'EH'+str(int(v_['15']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            vname = 'h'+str(int(v_['11']))+','+str(I)
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
            posev = np.where(elftv[ielftype[ie]]=='Y')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
        for I in range(int(v_['1']),int(v_['NSTEP'])+1):
            ename = 'PPW'+str(int(v_['11']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePMP4')
            ielftype = arrset(ielftype,ie,iet_["ePMP4"])
            ename = 'PPW'+str(int(v_['12']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePMP3')
            ielftype = arrset(ielftype,ie,iet_["ePMP3"])
            ename = 'PPW'+str(int(v_['13']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePMP1')
            ielftype = arrset(ielftype,ie,iet_["ePMP1"])
            ename = 'PPW'+str(int(v_['14']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePMP2')
            ielftype = arrset(ielftype,ie,iet_["ePMP2"])
            ename = 'PPW'+str(int(v_['15']))+','+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'ePMP2')
            ielftype = arrset(ielftype,ie,iet_["ePMP2"])
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                ename = 'PPW'+str(J)+','+str(I)
                [ie,ie_,_] = s2mpj_ii(ename,ie_)
                vname = 'q'+str(J)+','+str(I)
                [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,float(0.0))
                posev = np.where(elftv[ielftype[ie]]=='Q')[0]
                self.elvar = loaset(self.elvar,ie,posev[0],iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for I in range(int(v_['1']),int(v_['NSTEP'])+1):
            for J in range(int(v_['SuPIPE']),int(v_['EuPIPE'])+1):
                ig = ig_['EL'+str(J)+','+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['MGINV'+str(J)]))
            for J in range(int(v_['SuVAL']),int(v_['EuVAL'])+1):
                ig = ig_['EL'+str(J)+','+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['MGINV'+str(J)]))
            for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
                ig = ig_['EL'+str(J)+','+str(I)]
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['CONA'+str(J)]))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['CONB'+str(J)]))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(v_['CONC'+str(J)]))
                posel = len(self.grelt[ig])
                self.grelt = loaset(self.grelt,ig,posel,ie_['EH'+str(J)+','+str(I)])
                nlc = np.union1d(nlc,np.array([ig]))
                self.grelw = loaset(self.grelw,ig,posel,float(1.0))
        for J in range(int(v_['SuPMP']),int(v_['EuPMP'])+1):
            ig = ig_['PC'+str(J)+','+str(int(v_['2']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['PPW'+str(J)+','+str(int(v_['2']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(100.0))
            ig = ig_['PC'+str(J)+','+str(int(v_['3']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['PPW'+str(J)+','+str(int(v_['3']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(100.0))
            ig = ig_['PC'+str(J)+','+str(int(v_['4']))]
            posel = len(self.grelt[ig])
            self.grelt  = (
                  loaset(self.grelt,ig,posel,ie_['PPW'+str(J)+','+str(int(v_['4']))]))
            nlc = np.union1d(nlc,np.array([ig]))
            self.grelw = loaset(self.grelw,ig,posel,float(100.0))
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        self.A = csr_matrix((valA,(irA,icA)),shape=(ngrp,self.n))
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        self.clower = np.full((self.m,1),-float('Inf'))
        self.cupper = np.full((self.m,1),+float('Inf'))
        self.clower[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        self.cupper[np.arange(self.nle,self.nle+self.neq)] = np.zeros((self.neq,1))
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.lincons  = (
              np.where(np.isin(self.congrps,np.setdiff1d(self.congrps,nlc)))[0])
        self.pbclass   = "C-CSOR2-AY-312-256"
        self.objderlvl = 2
        self.conderlvl = [2]


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eXSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = EV_[0]*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = EV_[0]+EV_[0]
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXZ(self, nargout,*args):

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

    @staticmethod
    def eXMYZSQ(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        U_ = np.zeros((2,3))
        IV_ = np.zeros(2)
        U_[0,0] = U_[0,0]+1
        U_[0,1] = U_[0,1]-1
        U_[1,2] = U_[1,2]+1
        IV_[0] = U_[0:1,:].dot(EV_)
        IV_[1] = U_[1:2,:].dot(EV_)
        f_   = IV_[0]*IV_[1]**2
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = IV_[1]**2
            g_[1] = 2.0*IV_[0]*IV_[1]
            g_ =  U_.T.dot(g_)
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = 2.0*IV_[1]
                H_[1,0] = H_[0,1]
                H_[1,1] = 2.0*IV_[0]
                H_ = U_.T.dot(H_).dot(U_)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eXPOW(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        MODX = np.absolute(EV_[0])
        ISNEG = EV_[0]<0.0
        ISPOS = EV_[0]>=0.0
        if ISNEG!=0:
            SGN = -1.0
        if ISPOS!=0:
            SGN = +1.0
        f_   = SGN*MODX**1.852
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.852*MODX**0.852
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = SGN*1.577904*MODX**(-0.148)
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePMP1(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.074*EV_[0]*EV_[0]+3.062*EV_[0]+50.357
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.148*EV_[0]+3.062
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.148
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePMP2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.747*EV_[0]*EV_[0]-10.287*EV_[0]
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 1.494*EV_[0]-10.287
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 1.494
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePMP3(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.034*EV_[0]*EV_[0]+0.220*EV_[0]+6.685
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.068*EV_[0]+0.220
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.068
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def ePMP4(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        f_   = 0.079*EV_[0]*EV_[0]-2.761*EV_[0]+35.014
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = 0.158*EV_[0]-2.761
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = 0.158
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

