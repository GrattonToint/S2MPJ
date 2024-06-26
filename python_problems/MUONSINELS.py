from s2mpjlib import *
class  MUONSINELS(CUTEst_problem):

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MUONSINELS
#    *********
# 
#    ISIS Data fitting problem MUOSINE as a least-squares problem
# 
#    Fit: y = sin( b * x ) + e
# 
#    Source: fit to a sine using simplified muon data
#      from Mantid (http://www.mantidproject.org)
# 
#    SIF input: Nick Gould and Tyrone Rees, Dec 2015
# 
#    classification = "SUR2-MN-1-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = 'MUONSINELS'

    def __init__(self, *args): 
        import numpy as np
        nargin   = len(args)

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = {}
        ix_ = {}
        ig_ = {}
        v_['M'] = 512
        v_['N'] = 1
        v_['1'] = 1
        v_['ONE'] = 1.0
        v_['X1'] = 0.0
        v_['X2'] = 0.0122718
        v_['X3'] = 0.0245437
        v_['X4'] = 0.0368155
        v_['X5'] = 0.0490874
        v_['X6'] = 0.0613592
        v_['X7'] = 0.0736311
        v_['X8'] = 0.0859029
        v_['X9'] = 0.0981748
        v_['X10'] = 0.110447
        v_['X11'] = 0.122718
        v_['X12'] = 0.13499
        v_['X13'] = 0.147262
        v_['X14'] = 0.159534
        v_['X15'] = 0.171806
        v_['X16'] = 0.184078
        v_['X17'] = 0.19635
        v_['X18'] = 0.208621
        v_['X19'] = 0.220893
        v_['X20'] = 0.233165
        v_['X21'] = 0.245437
        v_['X22'] = 0.257709
        v_['X23'] = 0.269981
        v_['X24'] = 0.282252
        v_['X25'] = 0.294524
        v_['X26'] = 0.306796
        v_['X27'] = 0.319068
        v_['X28'] = 0.33134
        v_['X29'] = 0.343612
        v_['X30'] = 0.355884
        v_['X31'] = 0.368155
        v_['X32'] = 0.380427
        v_['X33'] = 0.392699
        v_['X34'] = 0.404971
        v_['X35'] = 0.417243
        v_['X36'] = 0.429515
        v_['X37'] = 0.441786
        v_['X38'] = 0.454058
        v_['X39'] = 0.46633
        v_['X40'] = 0.478602
        v_['X41'] = 0.490874
        v_['X42'] = 0.503146
        v_['X43'] = 0.515418
        v_['X44'] = 0.527689
        v_['X45'] = 0.539961
        v_['X46'] = 0.552233
        v_['X47'] = 0.564505
        v_['X48'] = 0.576777
        v_['X49'] = 0.589049
        v_['X50'] = 0.60132
        v_['X51'] = 0.613592
        v_['X52'] = 0.625864
        v_['X53'] = 0.638136
        v_['X54'] = 0.650408
        v_['X55'] = 0.66268
        v_['X56'] = 0.674952
        v_['X57'] = 0.687223
        v_['X58'] = 0.699495
        v_['X59'] = 0.711767
        v_['X60'] = 0.724039
        v_['X61'] = 0.736311
        v_['X62'] = 0.748583
        v_['X63'] = 0.760854
        v_['X64'] = 0.773126
        v_['X65'] = 0.785398
        v_['X66'] = 0.79767
        v_['X67'] = 0.809942
        v_['X68'] = 0.822214
        v_['X69'] = 0.834486
        v_['X70'] = 0.846757
        v_['X71'] = 0.859029
        v_['X72'] = 0.871301
        v_['X73'] = 0.883573
        v_['X74'] = 0.895845
        v_['X75'] = 0.908117
        v_['X76'] = 0.920388
        v_['X77'] = 0.93266
        v_['X78'] = 0.944932
        v_['X79'] = 0.957204
        v_['X80'] = 0.969476
        v_['X81'] = 0.981748
        v_['X82'] = 0.99402
        v_['X83'] = 1.00629
        v_['X84'] = 1.01856
        v_['X85'] = 1.03084
        v_['X86'] = 1.04311
        v_['X87'] = 1.05538
        v_['X88'] = 1.06765
        v_['X89'] = 1.07992
        v_['X90'] = 1.09219
        v_['X91'] = 1.10447
        v_['X92'] = 1.11674
        v_['X93'] = 1.12901
        v_['X94'] = 1.14128
        v_['X95'] = 1.15355
        v_['X96'] = 1.16583
        v_['X97'] = 1.1781
        v_['X98'] = 1.19037
        v_['X99'] = 1.20264
        v_['X100'] = 1.21491
        v_['X101'] = 1.22718
        v_['X102'] = 1.23946
        v_['X103'] = 1.25173
        v_['X104'] = 1.264
        v_['X105'] = 1.27627
        v_['X106'] = 1.28854
        v_['X107'] = 1.30082
        v_['X108'] = 1.31309
        v_['X109'] = 1.32536
        v_['X110'] = 1.33763
        v_['X111'] = 1.3499
        v_['X112'] = 1.36217
        v_['X113'] = 1.37445
        v_['X114'] = 1.38672
        v_['X115'] = 1.39899
        v_['X116'] = 1.41126
        v_['X117'] = 1.42353
        v_['X118'] = 1.43581
        v_['X119'] = 1.44808
        v_['X120'] = 1.46035
        v_['X121'] = 1.47262
        v_['X122'] = 1.48489
        v_['X123'] = 1.49717
        v_['X124'] = 1.50944
        v_['X125'] = 1.52171
        v_['X126'] = 1.53398
        v_['X127'] = 1.54625
        v_['X128'] = 1.55852
        v_['X129'] = 1.5708
        v_['X130'] = 1.58307
        v_['X131'] = 1.59534
        v_['X132'] = 1.60761
        v_['X133'] = 1.61988
        v_['X134'] = 1.63216
        v_['X135'] = 1.64443
        v_['X136'] = 1.6567
        v_['X137'] = 1.66897
        v_['X138'] = 1.68124
        v_['X139'] = 1.69351
        v_['X140'] = 1.70579
        v_['X141'] = 1.71806
        v_['X142'] = 1.73033
        v_['X143'] = 1.7426
        v_['X144'] = 1.75487
        v_['X145'] = 1.76715
        v_['X146'] = 1.77942
        v_['X147'] = 1.79169
        v_['X148'] = 1.80396
        v_['X149'] = 1.81623
        v_['X150'] = 1.82851
        v_['X151'] = 1.84078
        v_['X152'] = 1.85305
        v_['X153'] = 1.86532
        v_['X154'] = 1.87759
        v_['X155'] = 1.88986
        v_['X156'] = 1.90214
        v_['X157'] = 1.91441
        v_['X158'] = 1.92668
        v_['X159'] = 1.93895
        v_['X160'] = 1.95122
        v_['X161'] = 1.9635
        v_['X162'] = 1.97577
        v_['X163'] = 1.98804
        v_['X164'] = 2.00031
        v_['X165'] = 2.01258
        v_['X166'] = 2.02485
        v_['X167'] = 2.03713
        v_['X168'] = 2.0494
        v_['X169'] = 2.06167
        v_['X170'] = 2.07394
        v_['X171'] = 2.08621
        v_['X172'] = 2.09849
        v_['X173'] = 2.11076
        v_['X174'] = 2.12303
        v_['X175'] = 2.1353
        v_['X176'] = 2.14757
        v_['X177'] = 2.15984
        v_['X178'] = 2.17212
        v_['X179'] = 2.18439
        v_['X180'] = 2.19666
        v_['X181'] = 2.20893
        v_['X182'] = 2.2212
        v_['X183'] = 2.23348
        v_['X184'] = 2.24575
        v_['X185'] = 2.25802
        v_['X186'] = 2.27029
        v_['X187'] = 2.28256
        v_['X188'] = 2.29484
        v_['X189'] = 2.30711
        v_['X190'] = 2.31938
        v_['X191'] = 2.33165
        v_['X192'] = 2.34392
        v_['X193'] = 2.35619
        v_['X194'] = 2.36847
        v_['X195'] = 2.38074
        v_['X196'] = 2.39301
        v_['X197'] = 2.40528
        v_['X198'] = 2.41755
        v_['X199'] = 2.42983
        v_['X200'] = 2.4421
        v_['X201'] = 2.45437
        v_['X202'] = 2.46664
        v_['X203'] = 2.47891
        v_['X204'] = 2.49118
        v_['X205'] = 2.50346
        v_['X206'] = 2.51573
        v_['X207'] = 2.528
        v_['X208'] = 2.54027
        v_['X209'] = 2.55254
        v_['X210'] = 2.56482
        v_['X211'] = 2.57709
        v_['X212'] = 2.58936
        v_['X213'] = 2.60163
        v_['X214'] = 2.6139
        v_['X215'] = 2.62618
        v_['X216'] = 2.63845
        v_['X217'] = 2.65072
        v_['X218'] = 2.66299
        v_['X219'] = 2.67526
        v_['X220'] = 2.68753
        v_['X221'] = 2.69981
        v_['X222'] = 2.71208
        v_['X223'] = 2.72435
        v_['X224'] = 2.73662
        v_['X225'] = 2.74889
        v_['X226'] = 2.76117
        v_['X227'] = 2.77344
        v_['X228'] = 2.78571
        v_['X229'] = 2.79798
        v_['X230'] = 2.81025
        v_['X231'] = 2.82252
        v_['X232'] = 2.8348
        v_['X233'] = 2.84707
        v_['X234'] = 2.85934
        v_['X235'] = 2.87161
        v_['X236'] = 2.88388
        v_['X237'] = 2.89616
        v_['X238'] = 2.90843
        v_['X239'] = 2.9207
        v_['X240'] = 2.93297
        v_['X241'] = 2.94524
        v_['X242'] = 2.95751
        v_['X243'] = 2.96979
        v_['X244'] = 2.98206
        v_['X245'] = 2.99433
        v_['X246'] = 3.0066
        v_['X247'] = 3.01887
        v_['X248'] = 3.03115
        v_['X249'] = 3.04342
        v_['X250'] = 3.05569
        v_['X251'] = 3.06796
        v_['X252'] = 3.08023
        v_['X253'] = 3.09251
        v_['X254'] = 3.10478
        v_['X255'] = 3.11705
        v_['X256'] = 3.12932
        v_['X257'] = 3.14159
        v_['X258'] = 3.15386
        v_['X259'] = 3.16614
        v_['X260'] = 3.17841
        v_['X261'] = 3.19068
        v_['X262'] = 3.20295
        v_['X263'] = 3.21522
        v_['X264'] = 3.2275
        v_['X265'] = 3.23977
        v_['X266'] = 3.25204
        v_['X267'] = 3.26431
        v_['X268'] = 3.27658
        v_['X269'] = 3.28885
        v_['X270'] = 3.30113
        v_['X271'] = 3.3134
        v_['X272'] = 3.32567
        v_['X273'] = 3.33794
        v_['X274'] = 3.35021
        v_['X275'] = 3.36249
        v_['X276'] = 3.37476
        v_['X277'] = 3.38703
        v_['X278'] = 3.3993
        v_['X279'] = 3.41157
        v_['X280'] = 3.42385
        v_['X281'] = 3.43612
        v_['X282'] = 3.44839
        v_['X283'] = 3.46066
        v_['X284'] = 3.47293
        v_['X285'] = 3.4852
        v_['X286'] = 3.49748
        v_['X287'] = 3.50975
        v_['X288'] = 3.52202
        v_['X289'] = 3.53429
        v_['X290'] = 3.54656
        v_['X291'] = 3.55884
        v_['X292'] = 3.57111
        v_['X293'] = 3.58338
        v_['X294'] = 3.59565
        v_['X295'] = 3.60792
        v_['X296'] = 3.62019
        v_['X297'] = 3.63247
        v_['X298'] = 3.64474
        v_['X299'] = 3.65701
        v_['X300'] = 3.66928
        v_['X301'] = 3.68155
        v_['X302'] = 3.69383
        v_['X303'] = 3.7061
        v_['X304'] = 3.71837
        v_['X305'] = 3.73064
        v_['X306'] = 3.74291
        v_['X307'] = 3.75518
        v_['X308'] = 3.76746
        v_['X309'] = 3.77973
        v_['X310'] = 3.792
        v_['X311'] = 3.80427
        v_['X312'] = 3.81654
        v_['X313'] = 3.82882
        v_['X314'] = 3.84109
        v_['X315'] = 3.85336
        v_['X316'] = 3.86563
        v_['X317'] = 3.8779
        v_['X318'] = 3.89018
        v_['X319'] = 3.90245
        v_['X320'] = 3.91472
        v_['X321'] = 3.92699
        v_['X322'] = 3.93926
        v_['X323'] = 3.95153
        v_['X324'] = 3.96381
        v_['X325'] = 3.97608
        v_['X326'] = 3.98835
        v_['X327'] = 4.00062
        v_['X328'] = 4.01289
        v_['X329'] = 4.02517
        v_['X330'] = 4.03744
        v_['X331'] = 4.04971
        v_['X332'] = 4.06198
        v_['X333'] = 4.07425
        v_['X334'] = 4.08652
        v_['X335'] = 4.0988
        v_['X336'] = 4.11107
        v_['X337'] = 4.12334
        v_['X338'] = 4.13561
        v_['X339'] = 4.14788
        v_['X340'] = 4.16016
        v_['X341'] = 4.17243
        v_['X342'] = 4.1847
        v_['X343'] = 4.19697
        v_['X344'] = 4.20924
        v_['X345'] = 4.22152
        v_['X346'] = 4.23379
        v_['X347'] = 4.24606
        v_['X348'] = 4.25833
        v_['X349'] = 4.2706
        v_['X350'] = 4.28287
        v_['X351'] = 4.29515
        v_['X352'] = 4.30742
        v_['X353'] = 4.31969
        v_['X354'] = 4.33196
        v_['X355'] = 4.34423
        v_['X356'] = 4.35651
        v_['X357'] = 4.36878
        v_['X358'] = 4.38105
        v_['X359'] = 4.39332
        v_['X360'] = 4.40559
        v_['X361'] = 4.41786
        v_['X362'] = 4.43014
        v_['X363'] = 4.44241
        v_['X364'] = 4.45468
        v_['X365'] = 4.46695
        v_['X366'] = 4.47922
        v_['X367'] = 4.4915
        v_['X368'] = 4.50377
        v_['X369'] = 4.51604
        v_['X370'] = 4.52831
        v_['X371'] = 4.54058
        v_['X372'] = 4.55285
        v_['X373'] = 4.56513
        v_['X374'] = 4.5774
        v_['X375'] = 4.58967
        v_['X376'] = 4.60194
        v_['X377'] = 4.61421
        v_['X378'] = 4.62649
        v_['X379'] = 4.63876
        v_['X380'] = 4.65103
        v_['X381'] = 4.6633
        v_['X382'] = 4.67557
        v_['X383'] = 4.68785
        v_['X384'] = 4.70012
        v_['X385'] = 4.71239
        v_['X386'] = 4.72466
        v_['X387'] = 4.73693
        v_['X388'] = 4.7492
        v_['X389'] = 4.76148
        v_['X390'] = 4.77375
        v_['X391'] = 4.78602
        v_['X392'] = 4.79829
        v_['X393'] = 4.81056
        v_['X394'] = 4.82284
        v_['X395'] = 4.83511
        v_['X396'] = 4.84738
        v_['X397'] = 4.85965
        v_['X398'] = 4.87192
        v_['X399'] = 4.88419
        v_['X400'] = 4.89647
        v_['X401'] = 4.90874
        v_['X402'] = 4.92101
        v_['X403'] = 4.93328
        v_['X404'] = 4.94555
        v_['X405'] = 4.95783
        v_['X406'] = 4.9701
        v_['X407'] = 4.98237
        v_['X408'] = 4.99464
        v_['X409'] = 5.00691
        v_['X410'] = 5.01919
        v_['X411'] = 5.03146
        v_['X412'] = 5.04373
        v_['X413'] = 5.056
        v_['X414'] = 5.06827
        v_['X415'] = 5.08054
        v_['X416'] = 5.09282
        v_['X417'] = 5.10509
        v_['X418'] = 5.11736
        v_['X419'] = 5.12963
        v_['X420'] = 5.1419
        v_['X421'] = 5.15418
        v_['X422'] = 5.16645
        v_['X423'] = 5.17872
        v_['X424'] = 5.19099
        v_['X425'] = 5.20326
        v_['X426'] = 5.21553
        v_['X427'] = 5.22781
        v_['X428'] = 5.24008
        v_['X429'] = 5.25235
        v_['X430'] = 5.26462
        v_['X431'] = 5.27689
        v_['X432'] = 5.28917
        v_['X433'] = 5.30144
        v_['X434'] = 5.31371
        v_['X435'] = 5.32598
        v_['X436'] = 5.33825
        v_['X437'] = 5.35052
        v_['X438'] = 5.3628
        v_['X439'] = 5.37507
        v_['X440'] = 5.38734
        v_['X441'] = 5.39961
        v_['X442'] = 5.41188
        v_['X443'] = 5.42416
        v_['X444'] = 5.43643
        v_['X445'] = 5.4487
        v_['X446'] = 5.46097
        v_['X447'] = 5.47324
        v_['X448'] = 5.48552
        v_['X449'] = 5.49779
        v_['X450'] = 5.51006
        v_['X451'] = 5.52233
        v_['X452'] = 5.5346
        v_['X453'] = 5.54687
        v_['X454'] = 5.55915
        v_['X455'] = 5.57142
        v_['X456'] = 5.58369
        v_['X457'] = 5.59596
        v_['X458'] = 5.60823
        v_['X459'] = 5.62051
        v_['X460'] = 5.63278
        v_['X461'] = 5.64505
        v_['X462'] = 5.65732
        v_['X463'] = 5.66959
        v_['X464'] = 5.68186
        v_['X465'] = 5.69414
        v_['X466'] = 5.70641
        v_['X467'] = 5.71868
        v_['X468'] = 5.73095
        v_['X469'] = 5.74322
        v_['X470'] = 5.7555
        v_['X471'] = 5.76777
        v_['X472'] = 5.78004
        v_['X473'] = 5.79231
        v_['X474'] = 5.80458
        v_['X475'] = 5.81686
        v_['X476'] = 5.82913
        v_['X477'] = 5.8414
        v_['X478'] = 5.85367
        v_['X479'] = 5.86594
        v_['X480'] = 5.87821
        v_['X481'] = 5.89049
        v_['X482'] = 5.90276
        v_['X483'] = 5.91503
        v_['X484'] = 5.9273
        v_['X485'] = 5.93957
        v_['X486'] = 5.95185
        v_['X487'] = 5.96412
        v_['X488'] = 5.97639
        v_['X489'] = 5.98866
        v_['X490'] = 6.00093
        v_['X491'] = 6.0132
        v_['X492'] = 6.02548
        v_['X493'] = 6.03775
        v_['X494'] = 6.05002
        v_['X495'] = 6.06229
        v_['X496'] = 6.07456
        v_['X497'] = 6.08684
        v_['X498'] = 6.09911
        v_['X499'] = 6.11138
        v_['X500'] = 6.12365
        v_['X501'] = 6.13592
        v_['X502'] = 6.14819
        v_['X503'] = 6.16047
        v_['X504'] = 6.17274
        v_['X505'] = 6.18501
        v_['X506'] = 6.19728
        v_['X507'] = 6.20955
        v_['X508'] = 6.22183
        v_['X509'] = 6.2341
        v_['X510'] = 6.24637
        v_['X511'] = 6.25864
        v_['X512'] = 6.27091
        v_['Y1'] = 0.0
        v_['Y2'] = 0.0735646
        v_['Y3'] = 0.14673
        v_['Y4'] = 0.219101
        v_['Y5'] = 0.290285
        v_['Y6'] = 0.359895
        v_['Y7'] = 0.427555
        v_['Y8'] = 0.492898
        v_['Y9'] = 0.55557
        v_['Y10'] = 0.615232
        v_['Y11'] = 0.671559
        v_['Y12'] = 0.724247
        v_['Y13'] = 0.77301
        v_['Y14'] = 0.817585
        v_['Y15'] = 0.857729
        v_['Y16'] = 0.893224
        v_['Y17'] = 0.92388
        v_['Y18'] = 0.949528
        v_['Y19'] = 0.970031
        v_['Y20'] = 0.985278
        v_['Y21'] = 0.995185
        v_['Y22'] = 0.999699
        v_['Y23'] = 0.998795
        v_['Y24'] = 0.99248
        v_['Y25'] = 0.980785
        v_['Y26'] = 0.963776
        v_['Y27'] = 0.941544
        v_['Y28'] = 0.91421
        v_['Y29'] = 0.881921
        v_['Y30'] = 0.844854
        v_['Y31'] = 0.803208
        v_['Y32'] = 0.757209
        v_['Y33'] = 0.707107
        v_['Y34'] = 0.653173
        v_['Y35'] = 0.595699
        v_['Y36'] = 0.534998
        v_['Y37'] = 0.471397
        v_['Y38'] = 0.405241
        v_['Y39'] = 0.33689
        v_['Y40'] = 0.266713
        v_['Y41'] = 0.19509
        v_['Y42'] = 0.122411
        v_['Y43'] = 0.0490677
        v_['Y44'] = -0.0245412
        v_['Y45'] = -0.0980171
        v_['Y46'] = -0.170962
        v_['Y47'] = -0.24298
        v_['Y48'] = -0.313682
        v_['Y49'] = -0.382683
        v_['Y50'] = -0.449611
        v_['Y51'] = -0.514103
        v_['Y52'] = -0.575808
        v_['Y53'] = -0.634393
        v_['Y54'] = -0.689541
        v_['Y55'] = -0.740951
        v_['Y56'] = -0.788346
        v_['Y57'] = -0.83147
        v_['Y58'] = -0.870087
        v_['Y59'] = -0.903989
        v_['Y60'] = -0.932993
        v_['Y61'] = -0.95694
        v_['Y62'] = -0.975702
        v_['Y63'] = -0.989177
        v_['Y64'] = -0.99729
        v_['Y65'] = -1.0
        v_['Y66'] = -0.99729
        v_['Y67'] = -0.989177
        v_['Y68'] = -0.975702
        v_['Y69'] = -0.95694
        v_['Y70'] = -0.932993
        v_['Y71'] = -0.903989
        v_['Y72'] = -0.870087
        v_['Y73'] = -0.83147
        v_['Y74'] = -0.788346
        v_['Y75'] = -0.740951
        v_['Y76'] = -0.689541
        v_['Y77'] = -0.634393
        v_['Y78'] = -0.575808
        v_['Y79'] = -0.514103
        v_['Y80'] = -0.449611
        v_['Y81'] = -0.382683
        v_['Y82'] = -0.313682
        v_['Y83'] = -0.24298
        v_['Y84'] = -0.170962
        v_['Y85'] = -0.0980171
        v_['Y86'] = -0.0245412
        v_['Y87'] = 0.0490677
        v_['Y88'] = 0.122411
        v_['Y89'] = 0.19509
        v_['Y90'] = 0.266713
        v_['Y91'] = 0.33689
        v_['Y92'] = 0.405241
        v_['Y93'] = 0.471397
        v_['Y94'] = 0.534998
        v_['Y95'] = 0.595699
        v_['Y96'] = 0.653173
        v_['Y97'] = 0.707107
        v_['Y98'] = 0.757209
        v_['Y99'] = 0.803208
        v_['Y100'] = 0.844854
        v_['Y101'] = 0.881921
        v_['Y102'] = 0.91421
        v_['Y103'] = 0.941544
        v_['Y104'] = 0.963776
        v_['Y105'] = 0.980785
        v_['Y106'] = 0.99248
        v_['Y107'] = 0.998795
        v_['Y108'] = 0.999699
        v_['Y109'] = 0.995185
        v_['Y110'] = 0.985278
        v_['Y111'] = 0.970031
        v_['Y112'] = 0.949528
        v_['Y113'] = 0.92388
        v_['Y114'] = 0.893224
        v_['Y115'] = 0.857729
        v_['Y116'] = 0.817585
        v_['Y117'] = 0.77301
        v_['Y118'] = 0.724247
        v_['Y119'] = 0.671559
        v_['Y120'] = 0.615232
        v_['Y121'] = 0.55557
        v_['Y122'] = 0.492898
        v_['Y123'] = 0.427555
        v_['Y124'] = 0.359895
        v_['Y125'] = 0.290285
        v_['Y126'] = 0.219101
        v_['Y127'] = 0.14673
        v_['Y128'] = 0.0735646
        v_['Y129'] = 3.67394e-16
        v_['Y130'] = -0.0735646
        v_['Y131'] = -0.14673
        v_['Y132'] = -0.219101
        v_['Y133'] = -0.290285
        v_['Y134'] = -0.359895
        v_['Y135'] = -0.427555
        v_['Y136'] = -0.492898
        v_['Y137'] = -0.55557
        v_['Y138'] = -0.615232
        v_['Y139'] = -0.671559
        v_['Y140'] = -0.724247
        v_['Y141'] = -0.77301
        v_['Y142'] = -0.817585
        v_['Y143'] = -0.857729
        v_['Y144'] = -0.893224
        v_['Y145'] = -0.92388
        v_['Y146'] = -0.949528
        v_['Y147'] = -0.970031
        v_['Y148'] = -0.985278
        v_['Y149'] = -0.995185
        v_['Y150'] = -0.999699
        v_['Y151'] = -0.998795
        v_['Y152'] = -0.99248
        v_['Y153'] = -0.980785
        v_['Y154'] = -0.963776
        v_['Y155'] = -0.941544
        v_['Y156'] = -0.91421
        v_['Y157'] = -0.881921
        v_['Y158'] = -0.844854
        v_['Y159'] = -0.803208
        v_['Y160'] = -0.757209
        v_['Y161'] = -0.707107
        v_['Y162'] = -0.653173
        v_['Y163'] = -0.595699
        v_['Y164'] = -0.534998
        v_['Y165'] = -0.471397
        v_['Y166'] = -0.405241
        v_['Y167'] = -0.33689
        v_['Y168'] = -0.266713
        v_['Y169'] = -0.19509
        v_['Y170'] = -0.122411
        v_['Y171'] = -0.0490677
        v_['Y172'] = 0.0245412
        v_['Y173'] = 0.0980171
        v_['Y174'] = 0.170962
        v_['Y175'] = 0.24298
        v_['Y176'] = 0.313682
        v_['Y177'] = 0.382683
        v_['Y178'] = 0.449611
        v_['Y179'] = 0.514103
        v_['Y180'] = 0.575808
        v_['Y181'] = 0.634393
        v_['Y182'] = 0.689541
        v_['Y183'] = 0.740951
        v_['Y184'] = 0.788346
        v_['Y185'] = 0.83147
        v_['Y186'] = 0.870087
        v_['Y187'] = 0.903989
        v_['Y188'] = 0.932993
        v_['Y189'] = 0.95694
        v_['Y190'] = 0.975702
        v_['Y191'] = 0.989177
        v_['Y192'] = 0.99729
        v_['Y193'] = 1.0
        v_['Y194'] = 0.99729
        v_['Y195'] = 0.989177
        v_['Y196'] = 0.975702
        v_['Y197'] = 0.95694
        v_['Y198'] = 0.932993
        v_['Y199'] = 0.903989
        v_['Y200'] = 0.870087
        v_['Y201'] = 0.83147
        v_['Y202'] = 0.788346
        v_['Y203'] = 0.740951
        v_['Y204'] = 0.689541
        v_['Y205'] = 0.634393
        v_['Y206'] = 0.575808
        v_['Y207'] = 0.514103
        v_['Y208'] = 0.449611
        v_['Y209'] = 0.382683
        v_['Y210'] = 0.313682
        v_['Y211'] = 0.24298
        v_['Y212'] = 0.170962
        v_['Y213'] = 0.0980171
        v_['Y214'] = 0.0245412
        v_['Y215'] = -0.0490677
        v_['Y216'] = -0.122411
        v_['Y217'] = -0.19509
        v_['Y218'] = -0.266713
        v_['Y219'] = -0.33689
        v_['Y220'] = -0.405241
        v_['Y221'] = -0.471397
        v_['Y222'] = -0.534998
        v_['Y223'] = -0.595699
        v_['Y224'] = -0.653173
        v_['Y225'] = -0.707107
        v_['Y226'] = -0.757209
        v_['Y227'] = -0.803208
        v_['Y228'] = -0.844854
        v_['Y229'] = -0.881921
        v_['Y230'] = -0.91421
        v_['Y231'] = -0.941544
        v_['Y232'] = -0.963776
        v_['Y233'] = -0.980785
        v_['Y234'] = -0.99248
        v_['Y235'] = -0.998795
        v_['Y236'] = -0.999699
        v_['Y237'] = -0.995185
        v_['Y238'] = -0.985278
        v_['Y239'] = -0.970031
        v_['Y240'] = -0.949528
        v_['Y241'] = -0.92388
        v_['Y242'] = -0.893224
        v_['Y243'] = -0.857729
        v_['Y244'] = -0.817585
        v_['Y245'] = -0.77301
        v_['Y246'] = -0.724247
        v_['Y247'] = -0.671559
        v_['Y248'] = -0.615232
        v_['Y249'] = -0.55557
        v_['Y250'] = -0.492898
        v_['Y251'] = -0.427555
        v_['Y252'] = -0.359895
        v_['Y253'] = -0.290285
        v_['Y254'] = -0.219101
        v_['Y255'] = -0.14673
        v_['Y256'] = -0.0735646
        v_['Y257'] = -7.34788e-16
        v_['Y258'] = 0.0735646
        v_['Y259'] = 0.14673
        v_['Y260'] = 0.219101
        v_['Y261'] = 0.290285
        v_['Y262'] = 0.359895
        v_['Y263'] = 0.427555
        v_['Y264'] = 0.492898
        v_['Y265'] = 0.55557
        v_['Y266'] = 0.615232
        v_['Y267'] = 0.671559
        v_['Y268'] = 0.724247
        v_['Y269'] = 0.77301
        v_['Y270'] = 0.817585
        v_['Y271'] = 0.857729
        v_['Y272'] = 0.893224
        v_['Y273'] = 0.92388
        v_['Y274'] = 0.949528
        v_['Y275'] = 0.970031
        v_['Y276'] = 0.985278
        v_['Y277'] = 0.995185
        v_['Y278'] = 0.999699
        v_['Y279'] = 0.998795
        v_['Y280'] = 0.99248
        v_['Y281'] = 0.980785
        v_['Y282'] = 0.963776
        v_['Y283'] = 0.941544
        v_['Y284'] = 0.91421
        v_['Y285'] = 0.881921
        v_['Y286'] = 0.844854
        v_['Y287'] = 0.803208
        v_['Y288'] = 0.757209
        v_['Y289'] = 0.707107
        v_['Y290'] = 0.653173
        v_['Y291'] = 0.595699
        v_['Y292'] = 0.534998
        v_['Y293'] = 0.471397
        v_['Y294'] = 0.405241
        v_['Y295'] = 0.33689
        v_['Y296'] = 0.266713
        v_['Y297'] = 0.19509
        v_['Y298'] = 0.122411
        v_['Y299'] = 0.0490677
        v_['Y300'] = -0.0245412
        v_['Y301'] = -0.0980171
        v_['Y302'] = -0.170962
        v_['Y303'] = -0.24298
        v_['Y304'] = -0.313682
        v_['Y305'] = -0.382683
        v_['Y306'] = -0.449611
        v_['Y307'] = -0.514103
        v_['Y308'] = -0.575808
        v_['Y309'] = -0.634393
        v_['Y310'] = -0.689541
        v_['Y311'] = -0.740951
        v_['Y312'] = -0.788346
        v_['Y313'] = -0.83147
        v_['Y314'] = -0.870087
        v_['Y315'] = -0.903989
        v_['Y316'] = -0.932993
        v_['Y317'] = -0.95694
        v_['Y318'] = -0.975702
        v_['Y319'] = -0.989177
        v_['Y320'] = -0.99729
        v_['Y321'] = -1.0
        v_['Y322'] = -0.99729
        v_['Y323'] = -0.989177
        v_['Y324'] = -0.975702
        v_['Y325'] = -0.95694
        v_['Y326'] = -0.932993
        v_['Y327'] = -0.903989
        v_['Y328'] = -0.870087
        v_['Y329'] = -0.83147
        v_['Y330'] = -0.788346
        v_['Y331'] = -0.740951
        v_['Y332'] = -0.689541
        v_['Y333'] = -0.634393
        v_['Y334'] = -0.575808
        v_['Y335'] = -0.514103
        v_['Y336'] = -0.449611
        v_['Y337'] = -0.382683
        v_['Y338'] = -0.313682
        v_['Y339'] = -0.24298
        v_['Y340'] = -0.170962
        v_['Y341'] = -0.0980171
        v_['Y342'] = -0.0245412
        v_['Y343'] = 0.0490677
        v_['Y344'] = 0.122411
        v_['Y345'] = 0.19509
        v_['Y346'] = 0.266713
        v_['Y347'] = 0.33689
        v_['Y348'] = 0.405241
        v_['Y349'] = 0.471397
        v_['Y350'] = 0.534998
        v_['Y351'] = 0.595699
        v_['Y352'] = 0.653173
        v_['Y353'] = 0.707107
        v_['Y354'] = 0.757209
        v_['Y355'] = 0.803208
        v_['Y356'] = 0.844854
        v_['Y357'] = 0.881921
        v_['Y358'] = 0.91421
        v_['Y359'] = 0.941544
        v_['Y360'] = 0.963776
        v_['Y361'] = 0.980785
        v_['Y362'] = 0.99248
        v_['Y363'] = 0.998795
        v_['Y364'] = 0.999699
        v_['Y365'] = 0.995185
        v_['Y366'] = 0.985278
        v_['Y367'] = 0.970031
        v_['Y368'] = 0.949528
        v_['Y369'] = 0.92388
        v_['Y370'] = 0.893224
        v_['Y371'] = 0.857729
        v_['Y372'] = 0.817585
        v_['Y373'] = 0.77301
        v_['Y374'] = 0.724247
        v_['Y375'] = 0.671559
        v_['Y376'] = 0.615232
        v_['Y377'] = 0.55557
        v_['Y378'] = 0.492898
        v_['Y379'] = 0.427555
        v_['Y380'] = 0.359895
        v_['Y381'] = 0.290285
        v_['Y382'] = 0.219101
        v_['Y383'] = 0.14673
        v_['Y384'] = 0.0735646
        v_['Y385'] = 1.10218e-15
        v_['Y386'] = -0.0735646
        v_['Y387'] = -0.14673
        v_['Y388'] = -0.219101
        v_['Y389'] = -0.290285
        v_['Y390'] = -0.359895
        v_['Y391'] = -0.427555
        v_['Y392'] = -0.492898
        v_['Y393'] = -0.55557
        v_['Y394'] = -0.615232
        v_['Y395'] = -0.671559
        v_['Y396'] = -0.724247
        v_['Y397'] = -0.77301
        v_['Y398'] = -0.817585
        v_['Y399'] = -0.857729
        v_['Y400'] = -0.893224
        v_['Y401'] = -0.92388
        v_['Y402'] = -0.949528
        v_['Y403'] = -0.970031
        v_['Y404'] = -0.985278
        v_['Y405'] = -0.995185
        v_['Y406'] = -0.999699
        v_['Y407'] = -0.998795
        v_['Y408'] = -0.99248
        v_['Y409'] = -0.980785
        v_['Y410'] = -0.963776
        v_['Y411'] = -0.941544
        v_['Y412'] = -0.91421
        v_['Y413'] = -0.881921
        v_['Y414'] = -0.844854
        v_['Y415'] = -0.803208
        v_['Y416'] = -0.757209
        v_['Y417'] = -0.707107
        v_['Y418'] = -0.653173
        v_['Y419'] = -0.595699
        v_['Y420'] = -0.534998
        v_['Y421'] = -0.471397
        v_['Y422'] = -0.405241
        v_['Y423'] = -0.33689
        v_['Y424'] = -0.266713
        v_['Y425'] = -0.19509
        v_['Y426'] = -0.122411
        v_['Y427'] = -0.0490677
        v_['Y428'] = 0.0245412
        v_['Y429'] = 0.0980171
        v_['Y430'] = 0.170962
        v_['Y431'] = 0.24298
        v_['Y432'] = 0.313682
        v_['Y433'] = 0.382683
        v_['Y434'] = 0.449611
        v_['Y435'] = 0.514103
        v_['Y436'] = 0.575808
        v_['Y437'] = 0.634393
        v_['Y438'] = 0.689541
        v_['Y439'] = 0.740951
        v_['Y440'] = 0.788346
        v_['Y441'] = 0.83147
        v_['Y442'] = 0.870087
        v_['Y443'] = 0.903989
        v_['Y444'] = 0.932993
        v_['Y445'] = 0.95694
        v_['Y446'] = 0.975702
        v_['Y447'] = 0.989177
        v_['Y448'] = 0.99729
        v_['Y449'] = 1.0
        v_['Y450'] = 0.99729
        v_['Y451'] = 0.989177
        v_['Y452'] = 0.975702
        v_['Y453'] = 0.95694
        v_['Y454'] = 0.932993
        v_['Y455'] = 0.903989
        v_['Y456'] = 0.870087
        v_['Y457'] = 0.83147
        v_['Y458'] = 0.788346
        v_['Y459'] = 0.740951
        v_['Y460'] = 0.689541
        v_['Y461'] = 0.634393
        v_['Y462'] = 0.575808
        v_['Y463'] = 0.514103
        v_['Y464'] = 0.449611
        v_['Y465'] = 0.382683
        v_['Y466'] = 0.313682
        v_['Y467'] = 0.24298
        v_['Y468'] = 0.170962
        v_['Y469'] = 0.0980171
        v_['Y470'] = 0.0245412
        v_['Y471'] = -0.0490677
        v_['Y472'] = -0.122411
        v_['Y473'] = -0.19509
        v_['Y474'] = -0.266713
        v_['Y475'] = -0.33689
        v_['Y476'] = -0.405241
        v_['Y477'] = -0.471397
        v_['Y478'] = -0.534998
        v_['Y479'] = -0.595699
        v_['Y480'] = -0.653173
        v_['Y481'] = -0.707107
        v_['Y482'] = -0.757209
        v_['Y483'] = -0.803208
        v_['Y484'] = -0.844854
        v_['Y485'] = -0.881921
        v_['Y486'] = -0.91421
        v_['Y487'] = -0.941544
        v_['Y488'] = -0.963776
        v_['Y489'] = -0.980785
        v_['Y490'] = -0.99248
        v_['Y491'] = -0.998795
        v_['Y492'] = -0.999699
        v_['Y493'] = -0.995185
        v_['Y494'] = -0.985278
        v_['Y495'] = -0.970031
        v_['Y496'] = -0.949528
        v_['Y497'] = -0.92388
        v_['Y498'] = -0.893224
        v_['Y499'] = -0.857729
        v_['Y500'] = -0.817585
        v_['Y501'] = -0.77301
        v_['Y502'] = -0.724247
        v_['Y503'] = -0.671559
        v_['Y504'] = -0.615232
        v_['Y505'] = -0.55557
        v_['Y506'] = -0.492898
        v_['Y507'] = -0.427555
        v_['Y508'] = -0.359895
        v_['Y509'] = -0.290285
        v_['Y510'] = -0.219101
        v_['Y511'] = -0.14673
        v_['Y512'] = -0.0735646
        v_['E'] = 0.1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        self.xnames = np.array([])
        self.xscale = np.array([])
        intvars   = np.array([])
        binvars   = np.array([])
        [iv,ix_,_] = s2mpj_ii('B',ix_)
        self.xnames=arrset(self.xnames,iv,'B')
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
            v_['EINV'] = v_['ONE']/v_['E']
            v_['YOVERE'] = v_['EINV']*v_['Y'+str(I)]
            self.gconst = arrset(self.gconst,ig_['F'+str(I)],float(v_['YOVERE']))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        self.xlower = np.full((self.n,1),-float('Inf'))
        self.xupper = np.full((self.n,1),+float('Inf'))
        self.xlower = np.zeros((self.n,1))
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        self.x0 = np.zeros((self.n,1))
        self.x0[ix_['B']] = float(5.2)
        pass
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = {}
        elftv = []
        [it,iet_,_] = s2mpj_ii( 'eSINE', iet_)
        elftv = loaset(elftv,it,0,'B')
        elftp = []
        elftp = loaset(elftp,it,0,'X')
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = {}
        self.elftype = np.array([])
        ielftype     = np.array([])
        self.elvar   = []
        self.elpar   = []
        for I in range(int(v_['1']),int(v_['M'])+1):
            ename = 'G'+str(I)
            [ie,ie_,_] = s2mpj_ii(ename,ie_)
            self.elftype = arrset(self.elftype,ie,'eSINE')
            ielftype = arrset(ielftype, ie, iet_["eSINE"])
            vname = 'B'
            [iv,ix_] = s2mpj_nlx(self,vname,ix_,1,None,None,None)
            posev = np.where(elftv[ielftype[ie]]=='B')[0]
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
            v_['EINV'] = v_['ONE']/v_['E']
            ig = ig_['F'+str(I)]
            posel = len(self.grelt[ig])
            self.grelt = loaset(self.grelt,ig,posel,ie_['G'+str(I)])
            self.grelw = loaset(self.grelw,ig,posel,float(v_['EINV']))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        self.objlower = 0.0
#    Solution
# LO SOLTN
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        delattr( self, "A" )
        #%%%% RETURN VALUES FROM THE __INIT__ METHOD %%%%%%
        self.pbclass = "SUR2-MN-1-0"
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    @staticmethod
    def eSINE(self, nargout,*args):

        import numpy as np
        EV_  = args[0]
        iel_ = args[1]
        XB = self.elpar[iel_][0]*EV_[0]
        S = np.sin(XB)
        C = np.cos(XB)
        f_   = S
        if not isinstance( f_, float ):
            f_   = f_.item();
        if nargout>1:
            try:
                dim = len(IV_)
            except:
                dim = len(EV_)
            g_ = np.zeros(dim)
            g_[0] = self.elpar[iel_][0]*C
            if nargout>2:
                H_ = np.zeros((1,1))
                H_[0,0] = -self.elpar[iel_][0]*self.elpar[iel_][0]*S
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

