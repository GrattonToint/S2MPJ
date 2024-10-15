function ENSO(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ENSO
#    *********
# 
#    NIST Data fitting problem ENSO given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 + b2*cos( 2*pi*x/12 ) + b3*sin( 2*pi*x/12 ) 
#                      + b5*cos( 2*pi*x/b4 ) + b6*sin( 2*pi*x/b4 )
#                      + b8*cos( 2*pi*x/b7 ) + b9*sin( 2*pi*x/b7 ) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Kahaner, D., C. Moler, and S. Nash, (1989). 
#     Numerical Methods and Software.  
#     Englewood Cliffs, NJ: Prentice Hall, pp. 441-445.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-NOR2-MN-9-168"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ENSO"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 168
        v_["N"] = 9
        v_["1"] = 1
        v_["TWELVE"] = 12.0
        v_["PI/4"] = atan(1.0)
        v_["2PI"] = 8.0*v_["PI/4"]
        v_["2PIBY12"] = v_["2PI"]/v_["TWELVE"]
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
        v_["Y1"] = 12.9
        v_["Y2"] = 11.3
        v_["Y3"] = 10.6
        v_["Y4"] = 11.2
        v_["Y5"] = 10.9
        v_["Y6"] = 7.5
        v_["Y7"] = 7.7
        v_["Y8"] = 11.7
        v_["Y9"] = 12.9
        v_["Y10"] = 14.3
        v_["Y11"] = 10.9
        v_["Y12"] = 13.7
        v_["Y13"] = 17.1
        v_["Y14"] = 14.0
        v_["Y15"] = 15.3
        v_["Y16"] = 8.5
        v_["Y17"] = 5.7
        v_["Y18"] = 5.5
        v_["Y19"] = 7.6
        v_["Y20"] = 8.6
        v_["Y21"] = 7.3
        v_["Y22"] = 7.6
        v_["Y23"] = 12.7
        v_["Y24"] = 11.0
        v_["Y25"] = 12.7
        v_["Y26"] = 12.9
        v_["Y27"] = 13.0
        v_["Y28"] = 10.9
        v_["Y29"] = 10.4
        v_["Y30"] = 10.2
        v_["Y31"] = 8.0
        v_["Y32"] = 10.9
        v_["Y33"] = 13.6
        v_["Y34"] = 10.5
        v_["Y35"] = 9.2
        v_["Y36"] = 12.4
        v_["Y37"] = 12.7
        v_["Y38"] = 13.3
        v_["Y39"] = 10.1
        v_["Y40"] = 7.8
        v_["Y41"] = 4.8
        v_["Y42"] = 3.0
        v_["Y43"] = 2.5
        v_["Y44"] = 6.3
        v_["Y45"] = 9.7
        v_["Y46"] = 11.6
        v_["Y47"] = 8.6
        v_["Y48"] = 12.4
        v_["Y49"] = 10.5
        v_["Y50"] = 13.3
        v_["Y51"] = 10.4
        v_["Y52"] = 8.1
        v_["Y53"] = 3.7
        v_["Y54"] = 10.7
        v_["Y55"] = 5.1
        v_["Y56"] = 10.4
        v_["Y57"] = 10.9
        v_["Y58"] = 11.7
        v_["Y59"] = 11.4
        v_["Y60"] = 13.7
        v_["Y61"] = 14.1
        v_["Y62"] = 14.0
        v_["Y63"] = 12.5
        v_["Y64"] = 6.3
        v_["Y65"] = 9.6
        v_["Y66"] = 11.7
        v_["Y67"] = 5.0
        v_["Y68"] = 10.8
        v_["Y69"] = 12.7
        v_["Y70"] = 10.8
        v_["Y71"] = 11.8
        v_["Y72"] = 12.6
        v_["Y73"] = 15.7
        v_["Y74"] = 12.6
        v_["Y75"] = 14.8
        v_["Y76"] = 7.8
        v_["Y77"] = 7.1
        v_["Y78"] = 11.2
        v_["Y79"] = 8.1
        v_["Y80"] = 6.4
        v_["Y81"] = 5.2
        v_["Y82"] = 12.0
        v_["Y83"] = 10.2
        v_["Y84"] = 12.7
        v_["Y85"] = 10.2
        v_["Y86"] = 14.7
        v_["Y87"] = 12.2
        v_["Y88"] = 7.1
        v_["Y89"] = 5.7
        v_["Y90"] = 6.7
        v_["Y91"] = 3.9
        v_["Y92"] = 8.5
        v_["Y93"] = 8.3
        v_["Y94"] = 10.8
        v_["Y95"] = 16.7
        v_["Y96"] = 12.6
        v_["Y97"] = 12.5
        v_["Y98"] = 12.5
        v_["Y99"] = 9.8
        v_["Y100"] = 7.2
        v_["Y101"] = 4.1
        v_["Y102"] = 10.6
        v_["Y103"] = 10.1
        v_["Y104"] = 10.1
        v_["Y105"] = 11.9
        v_["Y106"] = 13.6
        v_["Y107"] = 16.3
        v_["Y108"] = 17.6
        v_["Y109"] = 15.5
        v_["Y110"] = 16.0
        v_["Y111"] = 15.2
        v_["Y112"] = 11.2
        v_["Y113"] = 14.3
        v_["Y114"] = 14.5
        v_["Y115"] = 8.5
        v_["Y116"] = 12.0
        v_["Y117"] = 12.7
        v_["Y118"] = 11.3
        v_["Y119"] = 14.5
        v_["Y120"] = 15.1
        v_["Y121"] = 10.4
        v_["Y122"] = 11.5
        v_["Y123"] = 13.4
        v_["Y124"] = 7.5
        v_["Y125"] = 0.6
        v_["Y126"] = 0.3
        v_["Y127"] = 5.5
        v_["Y128"] = 5.0
        v_["Y129"] = 4.6
        v_["Y130"] = 8.2
        v_["Y131"] = 9.9
        v_["Y132"] = 9.2
        v_["Y133"] = 12.5
        v_["Y134"] = 10.9
        v_["Y135"] = 9.9
        v_["Y136"] = 8.9
        v_["Y137"] = 7.6
        v_["Y138"] = 9.5
        v_["Y139"] = 8.4
        v_["Y140"] = 10.7
        v_["Y141"] = 13.6
        v_["Y142"] = 13.7
        v_["Y143"] = 13.7
        v_["Y144"] = 16.5
        v_["Y145"] = 16.8
        v_["Y146"] = 17.1
        v_["Y147"] = 15.4
        v_["Y148"] = 9.5
        v_["Y149"] = 6.1
        v_["Y150"] = 10.1
        v_["Y151"] = 9.3
        v_["Y152"] = 5.3
        v_["Y153"] = 11.2
        v_["Y154"] = 16.6
        v_["Y155"] = 15.6
        v_["Y156"] = 12.0
        v_["Y157"] = 11.5
        v_["Y158"] = 8.6
        v_["Y159"] = 13.8
        v_["Y160"] = 8.7
        v_["Y161"] = 8.6
        v_["Y162"] = 8.6
        v_["Y163"] = 8.7
        v_["Y164"] = 12.8
        v_["Y165"] = 13.2
        v_["Y166"] = 14.0
        v_["Y167"] = 13.4
        v_["Y168"] = 14.8
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
            v_["ARG"] = v_["2PIBY12"]*v_["X"*string(I)]
            v_["C"] = cos(v_["ARG"])
            v_["S"] = sin(v_["ARG"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
            iv = ix_["B1"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["B2"]
            pbm.A[ig,iv] += Float64(v_["C"])
            iv = ix_["B3"]
            pbm.A[ig,iv] += Float64(v_["S"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        legrps = findall(x->x=="<=",gtype)
        eqgrps = findall(x->x=="==",gtype)
        gegrps = findall(x->x==">=",gtype)
        pb.nle = length(legrps)
        pb.neq = length(eqgrps)
        pb.nge = length(gegrps)
        pb.m   = pb.nle+pb.neq+pb.nge
        pbm.congrps = [[legrps;eqgrps];gegrps]
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
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
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(11.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(11.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(0.5)
        end
        if haskey(ix_,"B4")
            pb.x0[ix_["B4"]] = Float64(40.0)
        else
            pb.y0[findfirst(x->x==ig_["B4"],pbm.congrps)] = Float64(40.0)
        end
        if haskey(ix_,"B5")
            pb.x0[ix_["B5"]] = Float64(-0.7)
        else
            pb.y0[findfirst(x->x==ig_["B5"],pbm.congrps)] = Float64(-0.7)
        end
        if haskey(ix_,"B6")
            pb.x0[ix_["B6"]] = Float64(-1.3)
        else
            pb.y0[findfirst(x->x==ig_["B6"],pbm.congrps)] = Float64(-1.3)
        end
        if haskey(ix_,"B7")
            pb.x0[ix_["B7"]] = Float64(25.0)
        else
            pb.y0[findfirst(x->x==ig_["B7"],pbm.congrps)] = Float64(25.0)
        end
        if haskey(ix_,"B8")
            pb.x0[ix_["B8"]] = Float64(-0.3)
        else
            pb.y0[findfirst(x->x==ig_["B8"],pbm.congrps)] = Float64(-0.3)
        end
        if haskey(ix_,"B9")
            pb.x0[ix_["B9"]] = Float64(1.4)
        else
            pb.y0[findfirst(x->x==ig_["B9"],pbm.congrps)] = Float64(1.4)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE8", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eE9", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE8")
            arrset(ielftype,ie,iet_["eE8"])
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "EB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE9")
            arrset(ielftype,ie,iet_["eE9"])
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "EC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE8")
            arrset(ielftype,ie,iet_["eE8"])
            vname = "B7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
            ename = "ED"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE9")
            arrset(ielftype,ie,iet_["eE9"])
            vname = "B7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B9"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ED"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-MN-9-168"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,8.0*atan(1.0e0))
        return pbm

    elseif action == "eE8"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V12 = EV_[1]*EV_[1]
        V13 = EV_[1]*V12
        V14 = V12*V12
        TPIX = pbm.efpar[1]*pbm.elpar[iel_][1]
        TPIXV1 = TPIX/EV_[1]
        C = cos(TPIXV1)
        S = sin(TPIXV1)
        TPIXS = TPIX*S
        f_   = EV_[2]*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = TPIXS*EV_[2]/V12
            g_[2] = C
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = (-pbm.efpar[1]*pbm.efpar[1]*EV_[2]*C*pbm.elpar[iel_][1]^2/
                     V14-2.0*TPIX*EV_[2]*S/V13)
                H_[1,2] = TPIXS/V12
                H_[2,1] = H_[1,2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eE9"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V12 = EV_[1]*EV_[1]
        V13 = EV_[1]*V12
        V14 = V12*V12
        TPIX = pbm.efpar[1]*pbm.elpar[iel_][1]
        TPIXV1 = TPIX/EV_[1]
        C = cos(TPIXV1)
        S = sin(TPIXV1)
        TPIXC = TPIX*C
        f_   = EV_[2]*S
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -TPIXC*EV_[2]/V12
            g_[2] = S
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      2.0*TPIX*EV_[2]*C/V13-pbm.efpar[1]*pbm.efpar[1]*EV_[2]*S*pbm.elpar[iel_][1]^2/V14)
                H_[1,2] = -TPIXC/V12
                H_[2,1] = H_[1,2]
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
            pbm.has_globs = [1,0]
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

