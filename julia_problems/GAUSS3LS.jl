function GAUSS3LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
#    classification = "C-SUR2-MN-8-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "GAUSS3LS"

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
        v_["Y1"] = 97.58776
        v_["Y2"] = 97.76344
        v_["Y3"] = 96.56705
        v_["Y4"] = 92.52037
        v_["Y5"] = 91.15097
        v_["Y6"] = 95.21728
        v_["Y7"] = 90.21355
        v_["Y8"] = 89.29235
        v_["Y9"] = 91.51479
        v_["Y10"] = 89.60965
        v_["Y11"] = 86.56187
        v_["Y12"] = 85.55315
        v_["Y13"] = 87.13053
        v_["Y14"] = 85.67938
        v_["Y15"] = 80.04849
        v_["Y16"] = 82.18922
        v_["Y17"] = 87.24078
        v_["Y18"] = 80.79401
        v_["Y19"] = 81.28564
        v_["Y20"] = 81.56932
        v_["Y21"] = 79.22703
        v_["Y22"] = 79.43259
        v_["Y23"] = 77.90174
        v_["Y24"] = 76.75438
        v_["Y25"] = 77.17338
        v_["Y26"] = 74.27296
        v_["Y27"] = 73.11830
        v_["Y28"] = 73.84732
        v_["Y29"] = 72.47746
        v_["Y30"] = 71.92128
        v_["Y31"] = 66.91962
        v_["Y32"] = 67.93554
        v_["Y33"] = 69.55841
        v_["Y34"] = 69.06592
        v_["Y35"] = 66.53371
        v_["Y36"] = 63.87094
        v_["Y37"] = 69.70526
        v_["Y38"] = 63.59295
        v_["Y39"] = 63.35509
        v_["Y40"] = 59.99747
        v_["Y41"] = 62.64843
        v_["Y42"] = 65.77345
        v_["Y43"] = 59.10141
        v_["Y44"] = 56.57750
        v_["Y45"] = 61.15313
        v_["Y46"] = 54.30767
        v_["Y47"] = 62.83535
        v_["Y48"] = 56.52957
        v_["Y49"] = 56.98427
        v_["Y50"] = 58.11459
        v_["Y51"] = 58.69576
        v_["Y52"] = 58.23322
        v_["Y53"] = 54.90490
        v_["Y54"] = 57.91442
        v_["Y55"] = 56.96629
        v_["Y56"] = 51.13831
        v_["Y57"] = 49.27123
        v_["Y58"] = 52.92668
        v_["Y59"] = 54.47693
        v_["Y60"] = 51.81710
        v_["Y61"] = 51.05401
        v_["Y62"] = 52.51731
        v_["Y63"] = 51.83710
        v_["Y64"] = 54.48196
        v_["Y65"] = 49.05859
        v_["Y66"] = 50.52315
        v_["Y67"] = 50.32755
        v_["Y68"] = 46.44419
        v_["Y69"] = 50.89281
        v_["Y70"] = 52.13203
        v_["Y71"] = 49.78741
        v_["Y72"] = 49.01637
        v_["Y73"] = 54.18198
        v_["Y74"] = 53.17456
        v_["Y75"] = 53.20827
        v_["Y76"] = 57.43459
        v_["Y77"] = 51.95282
        v_["Y78"] = 54.20282
        v_["Y79"] = 57.46687
        v_["Y80"] = 53.60268
        v_["Y81"] = 58.86728
        v_["Y82"] = 57.66652
        v_["Y83"] = 63.71034
        v_["Y84"] = 65.24244
        v_["Y85"] = 65.10878
        v_["Y86"] = 69.96313
        v_["Y87"] = 68.85475
        v_["Y88"] = 73.32574
        v_["Y89"] = 76.21241
        v_["Y90"] = 78.06311
        v_["Y91"] = 75.37701
        v_["Y92"] = 87.54449
        v_["Y93"] = 89.50588
        v_["Y94"] = 95.82098
        v_["Y95"] = 97.48390
        v_["Y96"] = 100.86070
        v_["Y97"] = 102.48510
        v_["Y98"] = 105.7311
        v_["Y99"] = 111.3489
        v_["Y100"] = 111.0305
        v_["Y101"] = 110.1920
        v_["Y102"] = 118.3581
        v_["Y103"] = 118.8086
        v_["Y104"] = 122.4249
        v_["Y105"] = 124.0953
        v_["Y106"] = 125.9337
        v_["Y107"] = 127.8533
        v_["Y108"] = 131.0361
        v_["Y109"] = 133.3343
        v_["Y110"] = 135.1278
        v_["Y111"] = 131.7113
        v_["Y112"] = 131.9151
        v_["Y113"] = 132.1107
        v_["Y114"] = 127.6898
        v_["Y115"] = 133.2148
        v_["Y116"] = 128.2296
        v_["Y117"] = 133.5902
        v_["Y118"] = 127.2539
        v_["Y119"] = 128.3482
        v_["Y120"] = 124.8694
        v_["Y121"] = 124.6031
        v_["Y122"] = 117.0648
        v_["Y123"] = 118.1966
        v_["Y124"] = 119.5408
        v_["Y125"] = 114.7946
        v_["Y126"] = 114.2780
        v_["Y127"] = 120.3484
        v_["Y128"] = 114.8647
        v_["Y129"] = 111.6514
        v_["Y130"] = 110.1826
        v_["Y131"] = 108.4461
        v_["Y132"] = 109.0571
        v_["Y133"] = 106.5308
        v_["Y134"] = 109.4691
        v_["Y135"] = 106.8709
        v_["Y136"] = 107.3192
        v_["Y137"] = 106.9000
        v_["Y138"] = 109.6526
        v_["Y139"] = 107.1602
        v_["Y140"] = 108.2509
        v_["Y141"] = 104.96310
        v_["Y142"] = 109.3601
        v_["Y143"] = 107.6696
        v_["Y144"] = 99.77286
        v_["Y145"] = 104.96440
        v_["Y146"] = 106.1376
        v_["Y147"] = 106.5816
        v_["Y148"] = 100.12860
        v_["Y149"] = 101.66910
        v_["Y150"] = 96.44254
        v_["Y151"] = 97.34169
        v_["Y152"] = 96.97412
        v_["Y153"] = 90.73460
        v_["Y154"] = 93.37949
        v_["Y155"] = 82.12331
        v_["Y156"] = 83.01657
        v_["Y157"] = 78.87360
        v_["Y158"] = 74.86971
        v_["Y159"] = 72.79341
        v_["Y160"] = 65.14744
        v_["Y161"] = 67.02127
        v_["Y162"] = 60.16136
        v_["Y163"] = 57.13996
        v_["Y164"] = 54.05769
        v_["Y165"] = 50.42265
        v_["Y166"] = 47.82430
        v_["Y167"] = 42.85748
        v_["Y168"] = 42.45495
        v_["Y169"] = 38.30808
        v_["Y170"] = 36.95794
        v_["Y171"] = 33.94543
        v_["Y172"] = 34.19017
        v_["Y173"] = 31.66097
        v_["Y174"] = 23.56172
        v_["Y175"] = 29.61143
        v_["Y176"] = 23.88765
        v_["Y177"] = 22.49812
        v_["Y178"] = 24.86901
        v_["Y179"] = 17.29481
        v_["Y180"] = 18.09291
        v_["Y181"] = 15.34813
        v_["Y182"] = 14.77997
        v_["Y183"] = 13.87832
        v_["Y184"] = 12.88891
        v_["Y185"] = 16.20763
        v_["Y186"] = 16.29024
        v_["Y187"] = 15.29712
        v_["Y188"] = 14.97839
        v_["Y189"] = 12.11330
        v_["Y190"] = 14.24168
        v_["Y191"] = 12.53824
        v_["Y192"] = 15.19818
        v_["Y193"] = 11.70478
        v_["Y194"] = 15.83745
        v_["Y195"] = 10.035850
        v_["Y196"] = 9.307574
        v_["Y197"] = 12.86800
        v_["Y198"] = 8.571671
        v_["Y199"] = 11.60415
        v_["Y200"] = 12.42772
        v_["Y201"] = 11.23627
        v_["Y202"] = 11.13198
        v_["Y203"] = 7.761117
        v_["Y204"] = 6.758250
        v_["Y205"] = 14.23375
        v_["Y206"] = 10.63876
        v_["Y207"] = 8.893581
        v_["Y208"] = 11.55398
        v_["Y209"] = 11.57221
        v_["Y210"] = 11.58347
        v_["Y211"] = 9.724857
        v_["Y212"] = 11.43854
        v_["Y213"] = 11.22636
        v_["Y214"] = 10.170150
        v_["Y215"] = 12.50765
        v_["Y216"] = 6.200494
        v_["Y217"] = 9.018902
        v_["Y218"] = 10.80557
        v_["Y219"] = 13.09591
        v_["Y220"] = 3.914033
        v_["Y221"] = 9.567723
        v_["Y222"] = 8.038338
        v_["Y223"] = 10.230960
        v_["Y224"] = 9.367358
        v_["Y225"] = 7.695937
        v_["Y226"] = 6.118552
        v_["Y227"] = 8.793192
        v_["Y228"] = 7.796682
        v_["Y229"] = 12.45064
        v_["Y230"] = 10.61601
        v_["Y231"] = 6.001000
        v_["Y232"] = 6.765096
        v_["Y233"] = 8.764652
        v_["Y234"] = 4.586417
        v_["Y235"] = 8.390782
        v_["Y236"] = 7.209201
        v_["Y237"] = 10.012090
        v_["Y238"] = 7.327461
        v_["Y239"] = 6.525136
        v_["Y240"] = 2.840065
        v_["Y241"] = 10.323710
        v_["Y242"] = 4.790035
        v_["Y243"] = 8.376431
        v_["Y244"] = 6.263980
        v_["Y245"] = 2.705892
        v_["Y246"] = 8.362109
        v_["Y247"] = 8.983507
        v_["Y248"] = 3.362469
        v_["Y249"] = 1.182678
        v_["Y250"] = 4.875312
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
        pb.x0[ix_["B1"]] = Float64(94.9)
        pb.x0[ix_["B2"]] = Float64(0.009)
        pb.x0[ix_["B3"]] = Float64(90.1)
        pb.x0[ix_["B4"]] = Float64(113.0)
        pb.x0[ix_["B5"]] = Float64(20.0)
        pb.x0[ix_["B6"]] = Float64(73.8)
        pb.x0[ix_["B7"]] = Float64(140.0)
        pb.x0[ix_["B8"]] = Float64(20.0)
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

