function GAUSS1LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GAUSS1LS
#    *********
# 
#    NIST Data fitting problem GAUSS1.
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
#    classification = "C-SUR2-MN-8-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "GAUSS1LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 250
        v_["N"] = 8
        v_["1"] = 1
        v_["X1"] = 1.0
        v_["X2"] = 2.0
        v_["X3"] = 3.0
        v_["X4"] = 4.0
        v_["X5"] = 5.0
        v_["X6"] = 6.0
        v_["X7"] = 7.0
        v_["X8"] = 8.0
        v_["X9"] = 9.0
        v_["X10"] = 10.0
        v_["X11"] = 11.0
        v_["X12"] = 12.0
        v_["X13"] = 13.0
        v_["X14"] = 14.0
        v_["X15"] = 15.0
        v_["X16"] = 16.0
        v_["X17"] = 17.0
        v_["X18"] = 18.0
        v_["X19"] = 19.0
        v_["X20"] = 20.0
        v_["X21"] = 21.0
        v_["X22"] = 22.0
        v_["X23"] = 23.0
        v_["X24"] = 24.0
        v_["X25"] = 25.0
        v_["X26"] = 26.0
        v_["X27"] = 27.0
        v_["X28"] = 28.0
        v_["X29"] = 29.0
        v_["X30"] = 30.0
        v_["X31"] = 31.0
        v_["X32"] = 32.0
        v_["X33"] = 33.0
        v_["X34"] = 34.0
        v_["X35"] = 35.0
        v_["X36"] = 36.0
        v_["X37"] = 37.0
        v_["X38"] = 38.0
        v_["X39"] = 39.0
        v_["X40"] = 40.0
        v_["X41"] = 41.0
        v_["X42"] = 42.0
        v_["X43"] = 43.0
        v_["X44"] = 44.0
        v_["X45"] = 45.0
        v_["X46"] = 46.0
        v_["X47"] = 47.0
        v_["X48"] = 48.0
        v_["X49"] = 49.0
        v_["X50"] = 50.0
        v_["X51"] = 51.0
        v_["X52"] = 52.0
        v_["X53"] = 53.0
        v_["X54"] = 54.0
        v_["X55"] = 55.0
        v_["X56"] = 56.0
        v_["X57"] = 57.0
        v_["X58"] = 58.0
        v_["X59"] = 59.0
        v_["X60"] = 60.0
        v_["X61"] = 61.0
        v_["X62"] = 62.0
        v_["X63"] = 63.0
        v_["X64"] = 64.0
        v_["X65"] = 65.0
        v_["X66"] = 66.0
        v_["X67"] = 67.0
        v_["X68"] = 68.0
        v_["X69"] = 69.0
        v_["X70"] = 70.0
        v_["X71"] = 71.0
        v_["X72"] = 72.0
        v_["X73"] = 73.0
        v_["X74"] = 74.0
        v_["X75"] = 75.0
        v_["X76"] = 76.0
        v_["X77"] = 77.0
        v_["X78"] = 78.0
        v_["X79"] = 79.0
        v_["X80"] = 80.0
        v_["X81"] = 81.0
        v_["X82"] = 82.0
        v_["X83"] = 83.0
        v_["X84"] = 84.0
        v_["X85"] = 85.0
        v_["X86"] = 86.0
        v_["X87"] = 87.0
        v_["X88"] = 88.0
        v_["X89"] = 89.0
        v_["X90"] = 90.0
        v_["X91"] = 91.0
        v_["X92"] = 92.0
        v_["X93"] = 93.0
        v_["X94"] = 94.0
        v_["X95"] = 95.0
        v_["X96"] = 96.0
        v_["X97"] = 97.0
        v_["X98"] = 98.0
        v_["X99"] = 99.0
        v_["X100"] = 100.0
        v_["X101"] = 101.0
        v_["X102"] = 102.0
        v_["X103"] = 103.0
        v_["X104"] = 104.0
        v_["X105"] = 105.0
        v_["X106"] = 106.0
        v_["X107"] = 107.0
        v_["X108"] = 108.0
        v_["X109"] = 109.0
        v_["X110"] = 110.0
        v_["X111"] = 111.0
        v_["X112"] = 112.0
        v_["X113"] = 113.0
        v_["X114"] = 114.0
        v_["X115"] = 115.0
        v_["X116"] = 116.0
        v_["X117"] = 117.0
        v_["X118"] = 118.0
        v_["X119"] = 119.0
        v_["X120"] = 120.0
        v_["X121"] = 121.0
        v_["X122"] = 122.0
        v_["X123"] = 123.0
        v_["X124"] = 124.0
        v_["X125"] = 125.0
        v_["X126"] = 126.0
        v_["X127"] = 127.0
        v_["X128"] = 128.0
        v_["X129"] = 129.0
        v_["X130"] = 130.0
        v_["X131"] = 131.0
        v_["X132"] = 132.0
        v_["X133"] = 133.0
        v_["X134"] = 134.0
        v_["X135"] = 135.0
        v_["X136"] = 136.0
        v_["X137"] = 137.0
        v_["X138"] = 138.0
        v_["X139"] = 139.0
        v_["X140"] = 140.0
        v_["X141"] = 141.0
        v_["X142"] = 142.0
        v_["X143"] = 143.0
        v_["X144"] = 144.0
        v_["X145"] = 145.0
        v_["X146"] = 146.0
        v_["X147"] = 147.0
        v_["X148"] = 148.0
        v_["X149"] = 149.0
        v_["X150"] = 150.0
        v_["X151"] = 151.0
        v_["X152"] = 152.0
        v_["X153"] = 153.0
        v_["X154"] = 154.0
        v_["X155"] = 155.0
        v_["X156"] = 156.0
        v_["X157"] = 157.0
        v_["X158"] = 158.0
        v_["X159"] = 159.0
        v_["X160"] = 160.0
        v_["X161"] = 161.0
        v_["X162"] = 162.0
        v_["X163"] = 163.0
        v_["X164"] = 164.0
        v_["X165"] = 165.0
        v_["X166"] = 166.0
        v_["X167"] = 167.0
        v_["X168"] = 168.0
        v_["X169"] = 169.0
        v_["X170"] = 170.0
        v_["X171"] = 171.0
        v_["X172"] = 172.0
        v_["X173"] = 173.0
        v_["X174"] = 174.0
        v_["X175"] = 175.0
        v_["X176"] = 176.0
        v_["X177"] = 177.0
        v_["X178"] = 178.0
        v_["X179"] = 179.0
        v_["X180"] = 180.0
        v_["X181"] = 181.0
        v_["X182"] = 182.0
        v_["X183"] = 183.0
        v_["X184"] = 184.0
        v_["X185"] = 185.0
        v_["X186"] = 186.0
        v_["X187"] = 187.0
        v_["X188"] = 188.0
        v_["X189"] = 189.0
        v_["X190"] = 190.0
        v_["X191"] = 191.0
        v_["X192"] = 192.0
        v_["X193"] = 193.0
        v_["X194"] = 194.0
        v_["X195"] = 195.0
        v_["X196"] = 196.0
        v_["X197"] = 197.0
        v_["X198"] = 198.0
        v_["X199"] = 199.0
        v_["X200"] = 200.0
        v_["X201"] = 201.0
        v_["X202"] = 202.0
        v_["X203"] = 203.0
        v_["X204"] = 204.0
        v_["X205"] = 205.0
        v_["X206"] = 206.0
        v_["X207"] = 207.0
        v_["X208"] = 208.0
        v_["X209"] = 209.0
        v_["X210"] = 210.0
        v_["X211"] = 211.0
        v_["X212"] = 212.0
        v_["X213"] = 213.0
        v_["X214"] = 214.0
        v_["X215"] = 215.0
        v_["X216"] = 216.0
        v_["X217"] = 217.0
        v_["X218"] = 218.0
        v_["X219"] = 219.0
        v_["X220"] = 220.0
        v_["X221"] = 221.0
        v_["X222"] = 222.0
        v_["X223"] = 223.0
        v_["X224"] = 224.0
        v_["X225"] = 225.0
        v_["X226"] = 226.0
        v_["X227"] = 227.0
        v_["X228"] = 228.0
        v_["X229"] = 229.0
        v_["X230"] = 230.0
        v_["X231"] = 231.0
        v_["X232"] = 232.0
        v_["X233"] = 233.0
        v_["X234"] = 234.0
        v_["X235"] = 235.0
        v_["X236"] = 236.0
        v_["X237"] = 237.0
        v_["X238"] = 238.0
        v_["X239"] = 239.0
        v_["X240"] = 240.0
        v_["X241"] = 241.0
        v_["X242"] = 242.0
        v_["X243"] = 243.0
        v_["X244"] = 244.0
        v_["X245"] = 245.0
        v_["X246"] = 246.0
        v_["X247"] = 247.0
        v_["X248"] = 248.0
        v_["X249"] = 249.0
        v_["X250"] = 250.0
        v_["Y1"] = 97.62227
        v_["Y2"] = 97.80724
        v_["Y3"] = 96.62247
        v_["Y4"] = 92.59022
        v_["Y5"] = 91.23869
        v_["Y6"] = 95.32704
        v_["Y7"] = 90.35040
        v_["Y8"] = 89.46235
        v_["Y9"] = 91.72520
        v_["Y10"] = 89.86916
        v_["Y11"] = 86.88076
        v_["Y12"] = 85.94360
        v_["Y13"] = 87.60686
        v_["Y14"] = 86.25839
        v_["Y15"] = 80.74976
        v_["Y16"] = 83.03551
        v_["Y17"] = 88.25837
        v_["Y18"] = 82.01316
        v_["Y19"] = 82.74098
        v_["Y20"] = 83.30034
        v_["Y21"] = 81.27850
        v_["Y22"] = 81.85506
        v_["Y23"] = 80.75195
        v_["Y24"] = 80.09573
        v_["Y25"] = 81.07633
        v_["Y26"] = 78.81542
        v_["Y27"] = 78.38596
        v_["Y28"] = 79.93386
        v_["Y29"] = 79.48474
        v_["Y30"] = 79.95942
        v_["Y31"] = 76.10691
        v_["Y32"] = 78.39830
        v_["Y33"] = 81.43060
        v_["Y34"] = 82.48867
        v_["Y35"] = 81.65462
        v_["Y36"] = 80.84323
        v_["Y37"] = 88.68663
        v_["Y38"] = 84.74438
        v_["Y39"] = 86.83934
        v_["Y40"] = 85.97739
        v_["Y41"] = 91.28509
        v_["Y42"] = 97.22411
        v_["Y43"] = 93.51733
        v_["Y44"] = 94.10159
        v_["Y45"] = 101.91760
        v_["Y46"] = 98.43134
        v_["Y47"] = 110.4214
        v_["Y48"] = 107.6628
        v_["Y49"] = 111.7288
        v_["Y50"] = 116.5115
        v_["Y51"] = 120.7609
        v_["Y52"] = 123.9553
        v_["Y53"] = 124.2437
        v_["Y54"] = 130.7996
        v_["Y55"] = 133.2960
        v_["Y56"] = 130.7788
        v_["Y57"] = 132.0565
        v_["Y58"] = 138.6584
        v_["Y59"] = 142.9252
        v_["Y60"] = 142.7215
        v_["Y61"] = 144.1249
        v_["Y62"] = 147.4377
        v_["Y63"] = 148.2647
        v_["Y64"] = 152.0519
        v_["Y65"] = 147.3863
        v_["Y66"] = 149.2074
        v_["Y67"] = 148.9537
        v_["Y68"] = 144.5876
        v_["Y69"] = 148.1226
        v_["Y70"] = 148.0144
        v_["Y71"] = 143.8893
        v_["Y72"] = 140.9088
        v_["Y73"] = 143.4434
        v_["Y74"] = 139.3938
        v_["Y75"] = 135.9878
        v_["Y76"] = 136.3927
        v_["Y77"] = 126.7262
        v_["Y78"] = 124.4487
        v_["Y79"] = 122.8647
        v_["Y80"] = 113.8557
        v_["Y81"] = 113.7037
        v_["Y82"] = 106.8407
        v_["Y83"] = 107.0034
        v_["Y84"] = 102.46290
        v_["Y85"] = 96.09296
        v_["Y86"] = 94.57555
        v_["Y87"] = 86.98824
        v_["Y88"] = 84.90154
        v_["Y89"] = 81.18023
        v_["Y90"] = 76.40117
        v_["Y91"] = 67.09200
        v_["Y92"] = 72.67155
        v_["Y93"] = 68.10848
        v_["Y94"] = 67.99088
        v_["Y95"] = 63.34094
        v_["Y96"] = 60.55253
        v_["Y97"] = 56.18687
        v_["Y98"] = 53.64482
        v_["Y99"] = 53.70307
        v_["Y100"] = 48.07893
        v_["Y101"] = 42.21258
        v_["Y102"] = 45.65181
        v_["Y103"] = 41.69728
        v_["Y104"] = 41.24946
        v_["Y105"] = 39.21349
        v_["Y106"] = 37.71696
        v_["Y107"] = 36.68395
        v_["Y108"] = 37.30393
        v_["Y109"] = 37.43277
        v_["Y110"] = 37.45012
        v_["Y111"] = 32.64648
        v_["Y112"] = 31.84347
        v_["Y113"] = 31.39951
        v_["Y114"] = 26.68912
        v_["Y115"] = 32.25323
        v_["Y116"] = 27.61008
        v_["Y117"] = 33.58649
        v_["Y118"] = 28.10714
        v_["Y119"] = 30.26428
        v_["Y120"] = 28.01648
        v_["Y121"] = 29.11021
        v_["Y122"] = 23.02099
        v_["Y123"] = 25.65091
        v_["Y124"] = 28.50295
        v_["Y125"] = 25.23701
        v_["Y126"] = 26.13828
        v_["Y127"] = 33.53260
        v_["Y128"] = 29.25195
        v_["Y129"] = 27.09847
        v_["Y130"] = 26.52999
        v_["Y131"] = 25.52401
        v_["Y132"] = 26.69218
        v_["Y133"] = 24.55269
        v_["Y134"] = 27.71763
        v_["Y135"] = 25.20297
        v_["Y136"] = 25.61483
        v_["Y137"] = 25.06893
        v_["Y138"] = 27.63930
        v_["Y139"] = 24.94851
        v_["Y140"] = 25.86806
        v_["Y141"] = 22.48183
        v_["Y142"] = 26.90045
        v_["Y143"] = 25.39919
        v_["Y144"] = 17.90614
        v_["Y145"] = 23.76039
        v_["Y146"] = 25.89689
        v_["Y147"] = 27.64231
        v_["Y148"] = 22.86101
        v_["Y149"] = 26.47003
        v_["Y150"] = 23.72888
        v_["Y151"] = 27.54334
        v_["Y152"] = 30.52683
        v_["Y153"] = 28.07261
        v_["Y154"] = 34.92815
        v_["Y155"] = 28.29194
        v_["Y156"] = 34.19161
        v_["Y157"] = 35.41207
        v_["Y158"] = 37.09336
        v_["Y159"] = 40.98330
        v_["Y160"] = 39.53923
        v_["Y161"] = 47.80123
        v_["Y162"] = 47.46305
        v_["Y163"] = 51.04166
        v_["Y164"] = 54.58065
        v_["Y165"] = 57.53001
        v_["Y166"] = 61.42089
        v_["Y167"] = 62.79032
        v_["Y168"] = 68.51455
        v_["Y169"] = 70.23053
        v_["Y170"] = 74.42776
        v_["Y171"] = 76.59911
        v_["Y172"] = 81.62053
        v_["Y173"] = 83.42208
        v_["Y174"] = 79.17451
        v_["Y175"] = 88.56985
        v_["Y176"] = 85.66525
        v_["Y177"] = 86.55502
        v_["Y178"] = 90.65907
        v_["Y179"] = 84.27290
        v_["Y180"] = 85.72220
        v_["Y181"] = 83.10702
        v_["Y182"] = 82.16884
        v_["Y183"] = 80.42568
        v_["Y184"] = 78.15692
        v_["Y185"] = 79.79691
        v_["Y186"] = 77.84378
        v_["Y187"] = 74.50327
        v_["Y188"] = 71.57289
        v_["Y189"] = 65.88031
        v_["Y190"] = 65.01385
        v_["Y191"] = 60.19582
        v_["Y192"] = 59.66726
        v_["Y193"] = 52.95478
        v_["Y194"] = 53.87792
        v_["Y195"] = 44.91274
        v_["Y196"] = 41.09909
        v_["Y197"] = 41.68018
        v_["Y198"] = 34.53379
        v_["Y199"] = 34.86419
        v_["Y200"] = 33.14787
        v_["Y201"] = 29.58864
        v_["Y202"] = 27.29462
        v_["Y203"] = 21.91439
        v_["Y204"] = 19.08159
        v_["Y205"] = 24.90290
        v_["Y206"] = 19.82341
        v_["Y207"] = 16.75551
        v_["Y208"] = 18.24558
        v_["Y209"] = 17.23549
        v_["Y210"] = 16.34934
        v_["Y211"] = 13.71285
        v_["Y212"] = 14.75676
        v_["Y213"] = 13.97169
        v_["Y214"] = 12.42867
        v_["Y215"] = 14.35519
        v_["Y216"] = 7.703309
        v_["Y217"] = 10.234410
        v_["Y218"] = 11.78315
        v_["Y219"] = 13.87768
        v_["Y220"] = 4.535700
        v_["Y221"] = 10.059280
        v_["Y222"] = 8.424824
        v_["Y223"] = 10.533120
        v_["Y224"] = 9.602255
        v_["Y225"] = 7.877514
        v_["Y226"] = 6.258121
        v_["Y227"] = 8.899865
        v_["Y228"] = 7.877754
        v_["Y229"] = 12.51191
        v_["Y230"] = 10.66205
        v_["Y231"] = 6.035400
        v_["Y232"] = 6.790655
        v_["Y233"] = 8.783535
        v_["Y234"] = 4.600288
        v_["Y235"] = 8.400915
        v_["Y236"] = 7.216561
        v_["Y237"] = 10.017410
        v_["Y238"] = 7.331278
        v_["Y239"] = 6.527863
        v_["Y240"] = 2.842001
        v_["Y241"] = 10.325070
        v_["Y242"] = 4.790995
        v_["Y243"] = 8.377101
        v_["Y244"] = 6.264445
        v_["Y245"] = 2.706213
        v_["Y246"] = 8.362329
        v_["Y247"] = 8.983658
        v_["Y248"] = 3.362571
        v_["Y249"] = 1.182746
        v_["Y250"] = 4.875359
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["B1"]] = Float64(97.0)
        pb.x0[ix_["B2"]] = Float64(0.009)
        pb.x0[ix_["B3"]] = Float64(100.0)
        pb.x0[ix_["B4"]] = Float64(65.0)
        pb.x0[ix_["B5"]] = Float64(20.0)
        pb.x0[ix_["B6"]] = Float64(70.0)
        pb.x0[ix_["B7"]] = Float64(178.0)
        pb.x0[ix_["B8"]] = Float64(16.5)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE2", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eE17", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE2")
            arrset(ielftype,ie,iet_["eE2"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "EB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE17")
            arrset(ielftype,ie,iet_["eE17"])
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "EC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE17")
            arrset(ielftype,ie,iet_["eE17"])
            vname = "B6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-MN-8-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eE2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E = exp(-EV_[2]*pbm.elpar[iel_][1])
        V1E = EV_[1]*E
        f_   = V1E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E
            g_[2] = -V1E*pbm.elpar[iel_][1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -pbm.elpar[iel_][1]*E
                H_[2,1] = H_[1,2]
                H_[2,2] = V1E*pbm.elpar[iel_][1]^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eE17"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V2MX = EV_[2]-pbm.elpar[iel_][1]
        V2MX2 = V2MX*V2MX
        TV2MX = 2.0*V2MX
        TV2MX2 = 2.0*V2MX2
        R = V2MX/EV_[3]
        A = -R*R
        E = exp(A)
        V32 = EV_[3]*EV_[3]
        V33 = EV_[3]*V32
        V1E = EV_[1]*E
        TV1E = 2.0*V1E
        TV2MXV = TV2MX/V32
        TV2MXW = TV2MX2/V32
        f_   = V1E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E
            g_[2] = -V1E*TV2MXV
            g_[3] = TV1E*V2MX2/V33
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = -E*TV2MXV
                H_[2,1] = H_[1,2]
                H_[1,3] = E*TV2MX2/V33
                H_[3,1] = H_[1,3]
                H_[2,2] = TV1E*(TV2MXW-1.0)/V32
                H_[2,3] = TV1E*TV2MX*(1.0-V2MX2/V32)/V33
                H_[3,2] = H_[2,3]
                H_[3,3] = TV1E*V2MX2*(TV2MXW-3.0)/EV_[3]^4
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",
                       "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",
                       "LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: action "*action*" unavailable for problem "*name*".jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

