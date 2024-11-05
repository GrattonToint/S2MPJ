from s2mpjlib import *
class  GAUSS3LS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GAUSS3LS
#    *********
# 
#    NIST Data fitting problem GAUSS3.
# 
#    Fit: y = b1*exp( -b2*x ) + b3*exp( -(x-b4)**2 / b5**2 )
#                             + b6*exp( -(x-b7)**2 / b8**2 ) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Rust, B., NIST (1996).
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CSUR2-MN-8-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Python by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'GAUSS3LS'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 250
        v_['N'] = 8
        v_['1'] = 1
        v_['X1'] = 1.0
        v_['X2'] = 2.0
        v_['X3'] = 3.0
        v_['X4'] = 4.0
        v_['X5'] = 5.0
        v_['X6'] = 6.0
        v_['X7'] = 7.0
        v_['X8'] = 8.0
        v_['X9'] = 9.0
        v_['X10'] = 10.0
        v_['X11'] = 11.0
        v_['X12'] = 12.0
        v_['X13'] = 13.0
        v_['X14'] = 14.0
        v_['X15'] = 15.0
        v_['X16'] = 16.0
        v_['X17'] = 17.0
        v_['X18'] = 18.0
        v_['X19'] = 19.0
        v_['X20'] = 20.0
        v_['X21'] = 21.0
        v_['X22'] = 22.0
        v_['X23'] = 23.0
        v_['X24'] = 24.0
        v_['X25'] = 25.0
        v_['X26'] = 26.0
        v_['X27'] = 27.0
        v_['X28'] = 28.0
        v_['X29'] = 29.0
        v_['X30'] = 30.0
        v_['X31'] = 31.0
        v_['X32'] = 32.0
        v_['X33'] = 33.0
        v_['X34'] = 34.0
        v_['X35'] = 35.0
        v_['X36'] = 36.0
        v_['X37'] = 37.0
        v_['X38'] = 38.0
        v_['X39'] = 39.0
        v_['X40'] = 40.0
        v_['X41'] = 41.0
        v_['X42'] = 42.0
        v_['X43'] = 43.0
        v_['X44'] = 44.0
        v_['X45'] = 45.0
        v_['X46'] = 46.0
        v_['X47'] = 47.0
        v_['X48'] = 48.0
        v_['X49'] = 49.0
        v_['X50'] = 50.0
        v_['X51'] = 51.0
        v_['X52'] = 52.0
        v_['X53'] = 53.0
        v_['X54'] = 54.0
        v_['X55'] = 55.0
        v_['X56'] = 56.0
        v_['X57'] = 57.0
        v_['X58'] = 58.0
        v_['X59'] = 59.0
        v_['X60'] = 60.0
        v_['X61'] = 61.0
        v_['X62'] = 62.0
        v_['X63'] = 63.0
        v_['X64'] = 64.0
        v_['X65'] = 65.0
        v_['X66'] = 66.0
        v_['X67'] = 67.0
        v_['X68'] = 68.0
        v_['X69'] = 69.0
        v_['X70'] = 70.0
        v_['X71'] = 71.0
        v_['X72'] = 72.0
        v_['X73'] = 73.0
        v_['X74'] = 74.0
        v_['X75'] = 75.0
        v_['X76'] = 76.0
        v_['X77'] = 77.0
        v_['X78'] = 78.0
        v_['X79'] = 79.0
        v_['X80'] = 80.0
        v_['X81'] = 81.0
        v_['X82'] = 82.0
        v_['X83'] = 83.0
        v_['X84'] = 84.0
        v_['X85'] = 85.0
        v_['X86'] = 86.0
        v_['X87'] = 87.0
        v_['X88'] = 88.0
        v_['X89'] = 89.0
        v_['X90'] = 90.0
        v_['X91'] = 91.0
        v_['X92'] = 92.0
        v_['X93'] = 93.0
        v_['X94'] = 94.0
        v_['X95'] = 95.0
        v_['X96'] = 96.0
        v_['X97'] = 97.0
        v_['X98'] = 98.0
        v_['X99'] = 99.0
        v_['X100'] = 100.0
        v_['X101'] = 101.0
        v_['X102'] = 102.0
        v_['X103'] = 103.0
        v_['X104'] = 104.0
        v_['X105'] = 105.0
        v_['X106'] = 106.0
        v_['X107'] = 107.0
        v_['X108'] = 108.0
        v_['X109'] = 109.0
        v_['X110'] = 110.0
        v_['X111'] = 111.0
        v_['X112'] = 112.0
        v_['X113'] = 113.0
        v_['X114'] = 114.0
        v_['X115'] = 115.0
        v_['X116'] = 116.0
        v_['X117'] = 117.0
        v_['X118'] = 118.0
        v_['X119'] = 119.0
        v_['X120'] = 120.0
        v_['X121'] = 121.0
        v_['X122'] = 122.0
        v_['X123'] = 123.0
        v_['X124'] = 124.0
        v_['X125'] = 125.0
        v_['X126'] = 126.0
        v_['X127'] = 127.0
        v_['X128'] = 128.0
        v_['X129'] = 129.0
        v_['X130'] = 130.0
        v_['X131'] = 131.0
        v_['X132'] = 132.0
        v_['X133'] = 133.0
        v_['X134'] = 134.0
        v_['X135'] = 135.0
        v_['X136'] = 136.0
        v_['X137'] = 137.0
        v_['X138'] = 138.0
        v_['X139'] = 139.0
        v_['X140'] = 140.0
        v_['X141'] = 141.0
        v_['X142'] = 142.0
        v_['X143'] = 143.0
        v_['X144'] = 144.0
        v_['X145'] = 145.0
        v_['X146'] = 146.0
        v_['X147'] = 147.0
        v_['X148'] = 148.0
        v_['X149'] = 149.0
        v_['X150'] = 150.0
        v_['X151'] = 151.0
        v_['X152'] = 152.0
        v_['X153'] = 153.0
        v_['X154'] = 154.0
        v_['X155'] = 155.0
        v_['X156'] = 156.0
        v_['X157'] = 157.0
        v_['X158'] = 158.0
        v_['X159'] = 159.0
        v_['X160'] = 160.0
        v_['X161'] = 161.0
        v_['X162'] = 162.0
        v_['X163'] = 163.0
        v_['X164'] = 164.0
        v_['X165'] = 165.0
        v_['X166'] = 166.0
        v_['X167'] = 167.0
        v_['X168'] = 168.0
        v_['X169'] = 169.0
        v_['X170'] = 170.0
        v_['X171'] = 171.0
        v_['X172'] = 172.0
        v_['X173'] = 173.0
        v_['X174'] = 174.0
        v_['X175'] = 175.0
        v_['X176'] = 176.0
        v_['X177'] = 177.0
        v_['X178'] = 178.0
        v_['X179'] = 179.0
        v_['X180'] = 180.0
        v_['X181'] = 181.0
        v_['X182'] = 182.0
        v_['X183'] = 183.0
        v_['X184'] = 184.0
        v_['X185'] = 185.0
        v_['X186'] = 186.0
        v_['X187'] = 187.0
        v_['X188'] = 188.0
        v_['X189'] = 189.0
        v_['X190'] = 190.0
        v_['X191'] = 191.0
        v_['X192'] = 192.0
        v_['X193'] = 193.0
        v_['X194'] = 194.0
        v_['X195'] = 195.0
        v_['X196'] = 196.0
        v_['X197'] = 197.0
        v_['X198'] = 198.0
        v_['X199'] = 199.0
        v_['X200'] = 200.0
        v_['X201'] = 201.0
        v_['X202'] = 202.0
        v_['X203'] = 203.0
        v_['X204'] = 204.0
        v_['X205'] = 205.0
        v_['X206'] = 206.0
        v_['X207'] = 207.0
        v_['X208'] = 208.0
        v_['X209'] = 209.0
        v_['X210'] = 210.0
        v_['X211'] = 211.0
        v_['X212'] = 212.0
        v_['X213'] = 213.0
        v_['X214'] = 214.0
        v_['X215'] = 215.0
        v_['X216'] = 216.0
        v_['X217'] = 217.0
        v_['X218'] = 218.0
        v_['X219'] = 219.0
        v_['X220'] = 220.0
        v_['X221'] = 221.0
        v_['X222'] = 222.0
        v_['X223'] = 223.0
        v_['X224'] = 224.0
        v_['X225'] = 225.0
        v_['X226'] = 226.0
        v_['X227'] = 227.0
        v_['X228'] = 228.0
        v_['X229'] = 229.0
        v_['X230'] = 230.0
        v_['X231'] = 231.0
        v_['X232'] = 232.0
        v_['X233'] = 233.0
        v_['X234'] = 234.0
        v_['X235'] = 235.0
        v_['X236'] = 236.0
        v_['X237'] = 237.0
        v_['X238'] = 238.0
        v_['X239'] = 239.0
        v_['X240'] = 240.0
        v_['X241'] = 241.0
        v_['X242'] = 242.0
        v_['X243'] = 243.0
        v_['X244'] = 244.0
        v_['X245'] = 245.0
        v_['X246'] = 246.0
        v_['X247'] = 247.0
        v_['X248'] = 248.0
        v_['X249'] = 249.0
        v_['X250'] = 250.0
        v_['Y1'] = 97.58776
        v_['Y2'] = 97.76344
        v_['Y3'] = 96.56705
        v_['Y4'] = 92.52037
        v_['Y5'] = 91.15097
        v_['Y6'] = 95.21728
        v_['Y7'] = 90.21355
        v_['Y8'] = 89.29235
        v_['Y9'] = 91.51479
        v_['Y10'] = 89.60965
        v_['Y11'] = 86.56187
        v_['Y12'] = 85.55315
        v_['Y13'] = 87.13053
        v_['Y14'] = 85.67938
        v_['Y15'] = 80.04849
        v_['Y16'] = 82.18922
        v_['Y17'] = 87.24078
        v_['Y18'] = 80.79401
        v_['Y19'] = 81.28564
        v_['Y20'] = 81.56932
        v_['Y21'] = 79.22703
        v_['Y22'] = 79.43259
        v_['Y23'] = 77.90174
        v_['Y24'] = 76.75438
        v_['Y25'] = 77.17338
        v_['Y26'] = 74.27296
        v_['Y27'] = 73.11830
        v_['Y28'] = 73.84732
        v_['Y29'] = 72.47746
        v_['Y30'] = 71.92128
        v_['Y31'] = 66.91962
        v_['Y32'] = 67.93554
        v_['Y33'] = 69.55841
        v_['Y34'] = 69.06592
        v_['Y35'] = 66.53371
        v_['Y36'] = 63.87094
        v_['Y37'] = 69.70526
        v_['Y38'] = 63.59295
        v_['Y39'] = 63.35509
        v_['Y40'] = 59.99747
        v_['Y41'] = 62.64843
        v_['Y42'] = 65.77345
        v_['Y43'] = 59.10141
        v_['Y44'] = 56.57750
        v_['Y45'] = 61.15313
        v_['Y46'] = 54.30767
        v_['Y47'] = 62.83535
        v_['Y48'] = 56.52957
        v_['Y49'] = 56.98427
        v_['Y50'] = 58.11459
        v_['Y51'] = 58.69576
        v_['Y52'] = 58.23322
        v_['Y53'] = 54.90490
        v_['Y54'] = 57.91442
        v_['Y55'] = 56.96629
        v_['Y56'] = 51.13831
        v_['Y57'] = 49.27123
        v_['Y58'] = 52.92668
        v_['Y59'] = 54.47693
        v_['Y60'] = 51.81710
        v_['Y61'] = 51.05401
        v_['Y62'] = 52.51731
        v_['Y63'] = 51.83710
        v_['Y64'] = 54.48196
        v_['Y65'] = 49.05859
        v_['Y66'] = 50.52315
        v_['Y67'] = 50.32755
        v_['Y68'] = 46.44419
        v_['Y69'] = 50.89281
        v_['Y70'] = 52.13203
        v_['Y71'] = 49.78741
        v_['Y72'] = 49.01637
        v_['Y73'] = 54.18198
        v_['Y74'] = 53.17456
        v_['Y75'] = 53.20827
        v_['Y76'] = 57.43459
        v_['Y77'] = 51.95282
        v_['Y78'] = 54.20282
        v_['Y79'] = 57.46687
        v_['Y80'] = 53.60268
        v_['Y81'] = 58.86728
        v_['Y82'] = 57.66652
        v_['Y83'] = 63.71034
        v_['Y84'] = 65.24244
        v_['Y85'] = 65.10878
        v_['Y86'] = 69.96313
        v_['Y87'] = 68.85475
        v_['Y88'] = 73.32574
        v_['Y89'] = 76.21241
        v_['Y90'] = 78.06311
        v_['Y91'] = 75.37701
        v_['Y92'] = 87.54449
        v_['Y93'] = 89.50588
        v_['Y94'] = 95.82098
        v_['Y95'] = 97.48390
        v_['Y96'] = 100.86070
        v_['Y97'] = 102.48510
        v_['Y98'] = 105.7311
        v_['Y99'] = 111.3489
        v_['Y100'] = 111.0305
        v_['Y101'] = 110.1920
        v_['Y102'] = 118.3581
        v_['Y103'] = 118.8086
        v_['Y104'] = 122.4249
        v_['Y105'] = 124.0953
        v_['Y106'] = 125.9337
        v_['Y107'] = 127.8533
        v_['Y108'] = 131.0361
        v_['Y109'] = 133.3343
        v_['Y110'] = 135.1278
        v_['Y111'] = 131.7113
        v_['Y112'] = 131.9151
        v_['Y113'] = 132.1107
        v_['Y114'] = 127.6898
        v_['Y115'] = 133.2148
        v_['Y116'] = 128.2296
        v_['Y117'] = 133.5902
        v_['Y118'] = 127.2539
        v_['Y119'] = 128.3482
        v_['Y120'] = 124.8694
        v_['Y121'] = 124.6031
        v_['Y122'] = 117.0648
        v_['Y123'] = 118.1966
        v_['Y124'] = 119.5408
        v_['Y125'] = 114.7946
        v_['Y126'] = 114.2780
        v_['Y127'] = 120.3484
        v_['Y128'] = 114.8647
        v_['Y129'] = 111.6514
        v_['Y130'] = 110.1826
        v_['Y131'] = 108.4461
        v_['Y132'] = 109.0571
        v_['Y133'] = 106.5308
        v_['Y134'] = 109.4691
        v_['Y135'] = 106.8709
        v_['Y136'] = 107.3192
        v_['Y137'] = 106.9000
        v_['Y138'] = 109.6526
        v_['Y139'] = 107.1602
        v_['Y140'] = 108.2509
        v_['Y141'] = 104.96310
        v_['Y142'] = 109.3601
        v_['Y143'] = 107.6696
        v_['Y144'] = 99.77286
        v_['Y145'] = 104.96440
        v_['Y146'] = 106.1376
        v_['Y147'] = 106.5816
        v_['Y148'] = 100.12860
        v_['Y149'] = 101.66910
        v_['Y150'] = 96.44254
        v_['Y151'] = 97.34169
        v_['Y152'] = 96.97412
        v_['Y153'] = 90.73460
        v_['Y154'] = 93.37949
        v_['Y155'] = 82.12331
        v_['Y156'] = 83.01657
        v_['Y157'] = 78.87360
        v_['Y158'] = 74.86971
        v_['Y159'] = 72.79341
        v_['Y160'] = 65.14744
        v_['Y161'] = 67.02127
        v_['Y162'] = 60.16136
        v_['Y163'] = 57.13996
        v_['Y164'] = 54.05769
        v_['Y165'] = 50.42265
        v_['Y166'] = 47.82430
        v_['Y167'] = 42.85748
        v_['Y168'] = 42.45495
        v_['Y169'] = 38.30808
        v_['Y170'] = 36.95794
        v_['Y171'] = 33.94543
        v_['Y172'] = 34.19017
        v_['Y173'] = 31.66097
        v_['Y174'] = 23.56172
        v_['Y175'] = 29.61143
        v_['Y176'] = 23.88765
        v_['Y177'] = 22.49812
        v_['Y178'] = 24.86901
        v_['Y179'] = 17.29481
        v_['Y180'] = 18.09291
        v_['Y181'] = 15.34813
        v_['Y182'] = 14.77997
        v_['Y183'] = 13.87832
        v_['Y184'] = 12.88891
        v_['Y185'] = 16.20763
        v_['Y186'] = 16.29024
        v_['Y187'] = 15.29712
        v_['Y188'] = 14.97839
        v_['Y189'] = 12.11330
        v_['Y190'] = 14.24168
        v_['Y191'] = 12.53824
        v_['Y192'] = 15.19818
        v_['Y193'] = 11.70478
        v_['Y194'] = 15.83745
        v_['Y195'] = 10.035850
        v_['Y196'] = 9.307574
        v_['Y197'] = 12.86800
        v_['Y198'] = 8.571671
        v_['Y199'] = 11.60415
        v_['Y200'] = 12.42772
        v_['Y201'] = 11.23627
        v_['Y202'] = 11.13198
        v_['Y203'] = 7.761117
        v_['Y204'] = 6.758250
        v_['Y205'] = 14.23375
        v_['Y206'] = 10.63876
        v_['Y207'] = 8.893581
        v_['Y208'] = 11.55398
        v_['Y209'] = 11.57221
        v_['Y210'] = 11.58347
        v_['Y211'] = 9.724857
        v_['Y212'] = 11.43854
        v_['Y213'] = 11.22636
        v_['Y214'] = 10.170150
        v_['Y215'] = 12.50765
        v_['Y216'] = 6.200494
        v_['Y217'] = 9.018902
        v_['Y218'] = 10.80557
        v_['Y219'] = 13.09591
        v_['Y220'] = 3.914033
        v_['Y221'] = 9.567723
        v_['Y222'] = 8.038338
        v_['Y223'] = 10.230960
        v_['Y224'] = 9.367358
        v_['Y225'] = 7.695937
        v_['Y226'] = 6.118552
        v_['Y227'] = 8.793192
        v_['Y228'] = 7.796682
        v_['Y229'] = 12.45064
        v_['Y230'] = 10.61601
        v_['Y231'] = 6.001000
        v_['Y232'] = 6.765096
        v_['Y233'] = 8.764652
        v_['Y234'] = 4.586417
        v_['Y235'] = 8.390782
        v_['Y236'] = 7.209201
        v_['Y237'] = 10.012090
        v_['Y238'] = 7.327461
        v_['Y239'] = 6.525136
        v_['Y240'] = 2.840065
        v_['Y241'] = 10.323710
        v_['Y242'] = 4.790035
        v_['Y243'] = 8.376431
        v_['Y244'] = 6.263980
        v_['Y245'] = 2.705892
        v_['Y246'] = 8.362109
        v_['Y247'] = 8.983507
        v_['Y248'] = 3.362469
        v_['Y249'] = 1.182678
        v_['Y250'] = 4.875312
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        for I in range(int(v_['1']),int(v_['N'])+1):
            [iv,ix_,_] = s2mpj_ii('B'+str(I),ix_)
            self.xnames=arrset(self.xnames,iv,'B'+str(I))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        self.A       = lil_matrix((1000000,1000000))
        self.gscale  = np.array([])
        self.grnames = np.array([])
        cnames      = np.array([])
        self.cnames = np.array([])
        gtype       = np.array([])
        for I in range(int(v_['1']),int(v_['M'])+1):
            [ig,ig_,_] = s2mpj_ii('F'+str(I),ig_)
            gtype = arrset(gtype,ig,'<>')
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        self.n   = len(ix_)
        ngrp   = len(ig_)
        self.objgrps = np.arange(ngrp)
        self.m       = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        self.gconst = np.zeros((ngrp,1))
        for I in range(int(v_['1']),int(v_['M'])+1):
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['Y'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['B1']] = float(94.9)
        self.x0[ix_['B2']] = float(0.009)
        self.x0[ix_['B3']] = float(90.1)
        self.x0[ix_['B4']] = float(113.0)
        self.x0[ix_['B5']] = float(20.0)
        self.x0[ix_['B6']] = float(73.8)
        self.x0[ix_['B7']] = float(140.0)
        self.x0[ix_['B8']] = float(20.0)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eE2', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        [it,iet_,_] = s2mpj_ii( 'eE17', iet_)
        elftv = loaset(elftv,it,0,'V1')
        elftv = loaset(elftv,it,1,'V2')
        elftv = loaset(elftv,it,2,'V3')
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'EA'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE2')
            ielftype = arrset(ielftype,ie,iet_["eE2"])
            vname = 'B1'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B2'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'EB'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE17')
            ielftype = arrset(ielftype,ie,iet_["eE17"])
            vname = 'B3'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B4'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B5'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
            ename = 'EC'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eE17')
            ielftype = arrset(ielftype,ie,iet_["eE17"])
            vname = 'B6'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V1')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B7'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V2')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            vname = 'B8'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='V3')[0]
            self.elvar = loaset(self.elvar,ie,posev[0],iv)
            posep = np.where(elftp[ielftype[ie]]=='X')[0]
            self.elpar = loaset(self.elpar,ie,posep[0],float(v_['X'+str(I)]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = {}
        [it,igt_,_] = s2mpj_ii('gL2',igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        self.grelt   = []
        for ig in np.arange(0,ngrp):
            self.grelt.append(np.array([]))
        self.grftype = np.array([])
        self.grelw   = []
        nlc         = np.array([])
        for ig in range(0,ngrp):
            self.grftype = arrset(self.grftype,ig,'gL2')
        for I in range(int(v_['1']),int(v_['M'])+1):
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EA'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EB'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['EC'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( self, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "C-CSUR2-MN-8-0"
        self.objderlvl = 2

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eE2(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        E = np.exp(-EV_[1]*self.elpar[iel_][0])
        V1E = EV_[0]*E
        f_   = V1E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E
            g_[1] = -V1E*self.elpar[iel_][0]
            if nargout>2:
                H_ = np.zeros((2,2))
                H_[0,1] = -self.elpar[iel_][0]*E
                H_[1,0] = H_[0,1]
                H_[1,1] = V1E*self.elpar[iel_][0]**2
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    @staticmethod
    def eE17(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        V2MX = EV_[1]-self.elpar[iel_][0]
        V2MX2 = V2MX*V2MX
        TV2MX = 2.0*V2MX
        TV2MX2 = 2.0*V2MX2
        R = V2MX/EV_[2]
        A = -R*R
        E = np.exp(A)
        V32 = EV_[2]*EV_[2]
        V33 = EV_[2]*V32
        V1E = EV_[0]*E
        TV1E = 2.0*V1E
        TV2MXV = TV2MX/V32
        TV2MXW = TV2MX2/V32
        f_   = V1E
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = E
            g_[1] = -V1E*TV2MXV
            g_[2] = TV1E*V2MX2/V33
            if nargout>2:
                H_ = np.zeros((3,3))
                H_[0,1] = -E*TV2MXV
                H_[1,0] = H_[0,1]
                H_[0,2] = E*TV2MX2/V33
                H_[2,0] = H_[0,2]
                H_[1,1] = TV1E*(TV2MXW-1.0)/V32
                H_[1,2] = TV1E*TV2MX*(1.0-V2MX2/V32)/V33
                H_[2,1] = H_[1,2]
                H_[2,2] = TV1E*V2MX2*(TV2MXW-3.0)/EV_[2]**4
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    @staticmethod
    def gL2(self,nargout,*args):

        GVAR_ = args[0]
        igr_  = args[1]
        f_= GVAR_*GVAR_
        if nargout>1:
            g_ = GVAR_+GVAR_
            if nargout>2:
                H_ = np.zeros((1,1))
                H_ = 2.0
        if nargout == 1:
            return f_
        elif nargout == 2:
            return f_,g_
        elif nargout == 3:
            return f_,g_,H_

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

