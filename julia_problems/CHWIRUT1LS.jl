function CHWIRUT1LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHWIRUT1LS
#    *********
# 
#    NIST Data fitting problem CHWIRUT1.
# 
#    Fit: y = exp[-b1*x]/(b2+b3*x) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Chwirut, D., NIST (197?).  
#      Ultrasonic Reference Block Study. 
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-SUR2-MN-3-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CHWIRUT1LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 214
        v_["N"] = 3
        v_["1"] = 1
        v_["X1"] = 0.5
        v_["X2"] = 0.625
        v_["X3"] = 0.75
        v_["X4"] = 0.875
        v_["X5"] = 1.0
        v_["X6"] = 1.25
        v_["X7"] = 1.75
        v_["X8"] = 2.25
        v_["X9"] = 1.75
        v_["X10"] = 2.25
        v_["X11"] = 2.75
        v_["X12"] = 3.25
        v_["X13"] = 3.75
        v_["X14"] = 4.25
        v_["X15"] = 4.75
        v_["X16"] = 5.25
        v_["X17"] = 5.75
        v_["X18"] = 0.5
        v_["X19"] = 0.625
        v_["X20"] = 0.75
        v_["X21"] = 0.875
        v_["X22"] = 1.0
        v_["X23"] = 1.25
        v_["X24"] = 1.75
        v_["X25"] = 2.25
        v_["X26"] = 1.75
        v_["X27"] = 2.25
        v_["X28"] = 2.75
        v_["X29"] = 3.25
        v_["X30"] = 3.75
        v_["X31"] = 4.25
        v_["X32"] = 4.75
        v_["X33"] = 5.25
        v_["X34"] = 5.75
        v_["X35"] = 0.5
        v_["X36"] = 0.625
        v_["X37"] = 0.75
        v_["X38"] = 0.875
        v_["X39"] = 1.0
        v_["X40"] = 1.25
        v_["X41"] = 1.75
        v_["X42"] = 2.25
        v_["X43"] = 1.75
        v_["X44"] = 2.25
        v_["X45"] = 2.75
        v_["X46"] = 3.25
        v_["X47"] = 3.75
        v_["X48"] = 4.25
        v_["X49"] = 4.75
        v_["X50"] = 5.25
        v_["X51"] = 5.75
        v_["X52"] = 0.5
        v_["X53"] = 0.625
        v_["X54"] = 0.75
        v_["X55"] = 0.875
        v_["X56"] = 1.0
        v_["X57"] = 1.25
        v_["X58"] = 1.75
        v_["X59"] = 2.25
        v_["X60"] = 1.75
        v_["X61"] = 2.25
        v_["X62"] = 2.75
        v_["X63"] = 3.25
        v_["X64"] = 3.75
        v_["X65"] = 4.25
        v_["X66"] = 4.75
        v_["X67"] = 5.25
        v_["X68"] = 5.75
        v_["X69"] = 0.5
        v_["X70"] = 0.75
        v_["X71"] = 1.5
        v_["X72"] = 3.0
        v_["X73"] = 3.0
        v_["X74"] = 3.0
        v_["X75"] = 6.0
        v_["X76"] = 0.5
        v_["X77"] = 0.75
        v_["X78"] = 1.5
        v_["X79"] = 3.0
        v_["X80"] = 3.0
        v_["X81"] = 3.0
        v_["X82"] = 6.0
        v_["X83"] = 0.5
        v_["X84"] = 0.75
        v_["X85"] = 1.5
        v_["X86"] = 3.0
        v_["X87"] = 3.0
        v_["X88"] = 3.0
        v_["X89"] = 6.0
        v_["X90"] = 0.5
        v_["X91"] = 0.75
        v_["X92"] = 1.5
        v_["X93"] = 3.0
        v_["X94"] = 6.0
        v_["X95"] = 3.0
        v_["X96"] = 3.0
        v_["X97"] = 6.0
        v_["X98"] = 0.5
        v_["X99"] = 0.75
        v_["X100"] = 1.0
        v_["X101"] = 1.5
        v_["X102"] = 2.0
        v_["X103"] = 2.0
        v_["X104"] = 2.5
        v_["X105"] = 3.0
        v_["X106"] = 4.0
        v_["X107"] = 5.0
        v_["X108"] = 6.0
        v_["X109"] = 0.5
        v_["X110"] = 0.75
        v_["X111"] = 1.0
        v_["X112"] = 1.5
        v_["X113"] = 2.0
        v_["X114"] = 2.0
        v_["X115"] = 2.5
        v_["X116"] = 3.0
        v_["X117"] = 4.0
        v_["X118"] = 5.0
        v_["X119"] = 6.0
        v_["X120"] = 0.5
        v_["X121"] = 0.75
        v_["X122"] = 1.0
        v_["X123"] = 1.5
        v_["X124"] = 2.0
        v_["X125"] = 2.0
        v_["X126"] = 2.5
        v_["X127"] = 3.0
        v_["X128"] = 4.0
        v_["X129"] = 5.0
        v_["X130"] = 6.0
        v_["X131"] = 0.5
        v_["X132"] = 0.625
        v_["X133"] = 0.75
        v_["X134"] = 0.875
        v_["X135"] = 1.0
        v_["X136"] = 1.25
        v_["X137"] = 2.25
        v_["X138"] = 2.25
        v_["X139"] = 2.75
        v_["X140"] = 3.25
        v_["X141"] = 3.75
        v_["X142"] = 4.25
        v_["X143"] = 4.75
        v_["X144"] = 5.25
        v_["X145"] = 5.75
        v_["X146"] = 3.0
        v_["X147"] = 3.0
        v_["X148"] = 3.0
        v_["X149"] = 3.0
        v_["X150"] = 3.0
        v_["X151"] = 3.0
        v_["X152"] = 0.5
        v_["X153"] = 0.75
        v_["X154"] = 1.0
        v_["X155"] = 1.5
        v_["X156"] = 2.0
        v_["X157"] = 2.5
        v_["X158"] = 2.0
        v_["X159"] = 2.5
        v_["X160"] = 3.0
        v_["X161"] = 4.0
        v_["X162"] = 5.0
        v_["X163"] = 6.0
        v_["X164"] = 0.5
        v_["X165"] = 0.75
        v_["X166"] = 1.0
        v_["X167"] = 1.5
        v_["X168"] = 2.0
        v_["X169"] = 2.5
        v_["X170"] = 2.0
        v_["X171"] = 2.5
        v_["X172"] = 3.0
        v_["X173"] = 4.0
        v_["X174"] = 5.0
        v_["X175"] = 6.0
        v_["X176"] = 0.5
        v_["X177"] = 0.75
        v_["X178"] = 1.0
        v_["X179"] = 1.5
        v_["X180"] = 2.0
        v_["X181"] = 2.5
        v_["X182"] = 2.0
        v_["X183"] = 2.5
        v_["X184"] = 3.0
        v_["X185"] = 4.0
        v_["X186"] = 5.0
        v_["X187"] = 6.0
        v_["X188"] = 3.0
        v_["X189"] = 0.5
        v_["X190"] = 0.75
        v_["X191"] = 1.5
        v_["X192"] = 3.0
        v_["X193"] = 6.0
        v_["X194"] = 3.0
        v_["X195"] = 6.0
        v_["X196"] = 3.0
        v_["X197"] = 3.0
        v_["X198"] = 3.0
        v_["X199"] = 1.75
        v_["X200"] = 1.75
        v_["X201"] = 0.5
        v_["X202"] = 0.75
        v_["X203"] = 1.75
        v_["X204"] = 1.75
        v_["X205"] = 2.75
        v_["X206"] = 3.75
        v_["X207"] = 1.75
        v_["X208"] = 1.75
        v_["X209"] = 0.5
        v_["X210"] = 0.75
        v_["X211"] = 2.75
        v_["X212"] = 3.75
        v_["X213"] = 1.75
        v_["X214"] = 1.75
        v_["Y1"] = 92.9
        v_["Y2"] = 78.7
        v_["Y3"] = 64.2
        v_["Y4"] = 64.9
        v_["Y5"] = 57.1
        v_["Y6"] = 43.3
        v_["Y7"] = 31.1
        v_["Y8"] = 23.6
        v_["Y9"] = 31.05
        v_["Y10"] = 23.7750
        v_["Y11"] = 17.7375
        v_["Y12"] = 13.8
        v_["Y13"] = 11.5875
        v_["Y14"] = 9.4125
        v_["Y15"] = 7.7250
        v_["Y16"] = 7.35
        v_["Y17"] = 8.0250
        v_["Y18"] = 90.6
        v_["Y19"] = 76.9
        v_["Y20"] = 71.6
        v_["Y21"] = 63.6
        v_["Y22"] = 54.0
        v_["Y23"] = 39.2
        v_["Y24"] = 29.3
        v_["Y25"] = 21.4
        v_["Y26"] = 29.1750
        v_["Y27"] = 22.1250
        v_["Y28"] = 17.5125
        v_["Y29"] = 14.25
        v_["Y30"] = 9.45
        v_["Y31"] = 9.15
        v_["Y32"] = 7.9125
        v_["Y33"] = 8.4750
        v_["Y34"] = 6.1125
        v_["Y35"] = 80.0
        v_["Y36"] = 79.0
        v_["Y37"] = 63.8
        v_["Y38"] = 57.2
        v_["Y39"] = 53.2
        v_["Y40"] = 42.5
        v_["Y41"] = 26.8
        v_["Y42"] = 20.4
        v_["Y43"] = 26.85
        v_["Y44"] = 21.0
        v_["Y45"] = 16.4625
        v_["Y46"] = 12.5250
        v_["Y47"] = 10.5375
        v_["Y48"] = 8.5875
        v_["Y49"] = 7.1250
        v_["Y50"] = 6.1125
        v_["Y51"] = 5.9625
        v_["Y52"] = 74.1
        v_["Y53"] = 67.3
        v_["Y54"] = 60.8
        v_["Y55"] = 55.5
        v_["Y56"] = 50.3
        v_["Y57"] = 41.0
        v_["Y58"] = 29.4
        v_["Y59"] = 20.4
        v_["Y60"] = 29.3625
        v_["Y61"] = 21.15
        v_["Y62"] = 16.7625
        v_["Y63"] = 13.2
        v_["Y64"] = 10.8750
        v_["Y65"] = 8.1750
        v_["Y66"] = 7.35
        v_["Y67"] = 5.9625
        v_["Y68"] = 5.6250
        v_["Y69"] = 81.5
        v_["Y70"] = 62.4
        v_["Y71"] = 32.5
        v_["Y72"] = 12.41
        v_["Y73"] = 13.12
        v_["Y74"] = 15.56
        v_["Y75"] = 5.63
        v_["Y76"] = 78.0
        v_["Y77"] = 59.9
        v_["Y78"] = 33.2
        v_["Y79"] = 13.84
        v_["Y80"] = 12.75
        v_["Y81"] = 14.62
        v_["Y82"] = 3.94
        v_["Y83"] = 76.8
        v_["Y84"] = 61.0
        v_["Y85"] = 32.9
        v_["Y86"] = 13.87
        v_["Y87"] = 11.81
        v_["Y88"] = 13.31
        v_["Y89"] = 5.44
        v_["Y90"] = 78.0
        v_["Y91"] = 63.5
        v_["Y92"] = 33.8
        v_["Y93"] = 12.56
        v_["Y94"] = 5.63
        v_["Y95"] = 12.75
        v_["Y96"] = 13.12
        v_["Y97"] = 5.44
        v_["Y98"] = 76.8
        v_["Y99"] = 60.0
        v_["Y100"] = 47.8
        v_["Y101"] = 32.0
        v_["Y102"] = 22.2
        v_["Y103"] = 22.57
        v_["Y104"] = 18.82
        v_["Y105"] = 13.95
        v_["Y106"] = 11.25
        v_["Y107"] = 9.0
        v_["Y108"] = 6.67
        v_["Y109"] = 75.8
        v_["Y110"] = 62.0
        v_["Y111"] = 48.8
        v_["Y112"] = 35.2
        v_["Y113"] = 20.0
        v_["Y114"] = 20.32
        v_["Y115"] = 19.31
        v_["Y116"] = 12.75
        v_["Y117"] = 10.42
        v_["Y118"] = 7.31
        v_["Y119"] = 7.42
        v_["Y120"] = 70.5
        v_["Y121"] = 59.5
        v_["Y122"] = 48.5
        v_["Y123"] = 35.8
        v_["Y124"] = 21.0
        v_["Y125"] = 21.67
        v_["Y126"] = 21.0
        v_["Y127"] = 15.64
        v_["Y128"] = 8.17
        v_["Y129"] = 8.55
        v_["Y130"] = 10.12
        v_["Y131"] = 78.0
        v_["Y132"] = 66.0
        v_["Y133"] = 62.0
        v_["Y134"] = 58.0
        v_["Y135"] = 47.7
        v_["Y136"] = 37.8
        v_["Y137"] = 20.2
        v_["Y138"] = 21.07
        v_["Y139"] = 13.87
        v_["Y140"] = 9.67
        v_["Y141"] = 7.76
        v_["Y142"] = 5.44
        v_["Y143"] = 4.87
        v_["Y144"] = 4.01
        v_["Y145"] = 3.75
        v_["Y146"] = 24.19
        v_["Y147"] = 25.76
        v_["Y148"] = 18.07
        v_["Y149"] = 11.81
        v_["Y150"] = 12.07
        v_["Y151"] = 16.12
        v_["Y152"] = 70.8
        v_["Y153"] = 54.7
        v_["Y154"] = 48.0
        v_["Y155"] = 39.8
        v_["Y156"] = 29.8
        v_["Y157"] = 23.7
        v_["Y158"] = 29.62
        v_["Y159"] = 23.81
        v_["Y160"] = 17.7
        v_["Y161"] = 11.55
        v_["Y162"] = 12.07
        v_["Y163"] = 8.74
        v_["Y164"] = 80.7
        v_["Y165"] = 61.3
        v_["Y166"] = 47.5
        v_["Y167"] = 29.0
        v_["Y168"] = 24.0
        v_["Y169"] = 17.7
        v_["Y170"] = 24.56
        v_["Y171"] = 18.67
        v_["Y172"] = 16.24
        v_["Y173"] = 8.74
        v_["Y174"] = 7.87
        v_["Y175"] = 8.51
        v_["Y176"] = 66.7
        v_["Y177"] = 59.2
        v_["Y178"] = 40.8
        v_["Y179"] = 30.7
        v_["Y180"] = 25.7
        v_["Y181"] = 16.3
        v_["Y182"] = 25.99
        v_["Y183"] = 16.95
        v_["Y184"] = 13.35
        v_["Y185"] = 8.62
        v_["Y186"] = 7.2
        v_["Y187"] = 6.64
        v_["Y188"] = 13.69
        v_["Y189"] = 81.0
        v_["Y190"] = 64.5
        v_["Y191"] = 35.5
        v_["Y192"] = 13.31
        v_["Y193"] = 4.87
        v_["Y194"] = 12.94
        v_["Y195"] = 5.06
        v_["Y196"] = 15.19
        v_["Y197"] = 14.62
        v_["Y198"] = 15.64
        v_["Y199"] = 25.5
        v_["Y200"] = 25.95
        v_["Y201"] = 81.7
        v_["Y202"] = 61.6
        v_["Y203"] = 29.8
        v_["Y204"] = 29.81
        v_["Y205"] = 17.17
        v_["Y206"] = 10.39
        v_["Y207"] = 28.4
        v_["Y208"] = 28.69
        v_["Y209"] = 81.3
        v_["Y210"] = 60.9
        v_["Y211"] = 16.65
        v_["Y212"] = 10.05
        v_["Y213"] = 28.9
        v_["Y214"] = 28.95
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
        pb.x0[ix_["B1"]] = Float64(0.1)
        pb.x0[ix_["B2"]] = Float64(0.01)
        pb.x0[ix_["B3"]] = Float64(0.02)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE16", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE16")
            arrset(ielftype,ie,iet_["eE16"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B3"
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
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
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
        pb.pbclass = "C-SUR2-MN-3-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eE16"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E = exp(-EV_[1]*pbm.elpar[iel_][1])
        EX = E*pbm.elpar[iel_][1]
        EX2 = EX*pbm.elpar[iel_][1]
        V2PV3X = EV_[2]+EV_[3]*pbm.elpar[iel_][1]
        V2PV32 = V2PV3X*V2PV3X
        V2PV33 = V2PV3X*V2PV32
        f_   = E/V2PV3X
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EX/V2PV3X
            g_[2] = -E/V2PV32
            g_[3] = -EX/V2PV32
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = EX2/V2PV3X
                H_[1,2] = EX/V2PV32
                H_[2,1] = H_[1,2]
                H_[1,3] = EX2/V2PV32
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*E/V2PV33
                H_[2,3] = 2.0*EX/V2PV33
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*EX2/V2PV33
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

