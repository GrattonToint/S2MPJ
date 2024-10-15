function GROWTH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : GROWTH
#    *********
#    GROWTH problem in 3 variables
# 
#    Fit the observed growth g(n) from Gaussian Elimination
#    with complete pivoting to a function of the form
#         U1 * n ** ( U2 + LOG(n) * U3 )
# 
#    SIF input: Nick Gould, Nov, 1991.
# 
#    classification = "C-NOR2-AN-3-12"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "GROWTH"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 3
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("U1",ix_)
        arrset(pb.xnames,iv,"U1")
        iv,ix_,_ = s2mpj_ii("U2",ix_)
        arrset(pb.xnames,iv,"U2")
        iv,ix_,_ = s2mpj_ii("U3",ix_)
        arrset(pb.xnames,iv,"U3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("G8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G8")
        ig,ig_,_ = s2mpj_ii("G9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G9")
        ig,ig_,_ = s2mpj_ii("G10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G10")
        ig,ig_,_ = s2mpj_ii("G11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G11")
        ig,ig_,_ = s2mpj_ii("G12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G12")
        ig,ig_,_ = s2mpj_ii("G13",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G13")
        ig,ig_,_ = s2mpj_ii("G14",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G14")
        ig,ig_,_ = s2mpj_ii("G15",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G15")
        ig,ig_,_ = s2mpj_ii("G16",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G16")
        ig,ig_,_ = s2mpj_ii("G18",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G18")
        ig,ig_,_ = s2mpj_ii("G20",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G20")
        ig,ig_,_ = s2mpj_ii("G25",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G25")
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
        pbm.gconst[ig_["G8"]] = Float64(8.0)
        pbm.gconst[ig_["G9"]] = Float64(8.4305)
        pbm.gconst[ig_["G10"]] = Float64(9.5294)
        pbm.gconst[ig_["G11"]] = Float64(10.4627)
        pbm.gconst[ig_["G12"]] = Float64(12.0)
        pbm.gconst[ig_["G13"]] = Float64(13.0205)
        pbm.gconst[ig_["G14"]] = Float64(14.5949)
        pbm.gconst[ig_["G15"]] = Float64(16.1078)
        pbm.gconst[ig_["G16"]] = Float64(18.0596)
        pbm.gconst[ig_["G18"]] = Float64(20.4569)
        pbm.gconst[ig_["G20"]] = Float64(24.25)
        pbm.gconst[ig_["G25"]] = Float64(32.9863)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["U1"]] = Float64(100.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eFIT", iet_)
        loaset(elftv,it,1,"U1")
        loaset(elftv,it,2,"U2")
        loaset(elftv,it,3,"U3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"RN")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "G8"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(8.0))
        ename = "G9"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.0))
        ename = "G10"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(10.0))
        ename = "G11"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(11.0))
        ename = "G12"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(12.0))
        ename = "G13"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(13.0))
        ename = "G14"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(14.0))
        ename = "G15"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(15.0))
        ename = "G16"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(16.0))
        ename = "G18"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(18.0))
        ename = "G20"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(20.0))
        ename = "G25"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eFIT")
            arrset(ielftype,ie,iet_["eFIT"])
        end
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(25.0))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["G8"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G9"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G9"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G10"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G10"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G11"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G12"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G12"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G13"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G14"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G14"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G15"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G15"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G16"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G16"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G18"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G18"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G20"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G20"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["G25"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["G25"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN                0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-AN-3-12"
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

    elseif action == "eFIT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LOGRN = log(pbm.elpar[iel_][1])
        POWER = pbm.elpar[iel_][1]^(EV_[2]+LOGRN*EV_[3])
        f_   = EV_[1]*POWER
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = POWER
            g_[2] = EV_[1]*POWER*LOGRN
            g_[3] = EV_[1]*POWER*LOGRN^2
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 0.0
                H_[1,2] = POWER*LOGRN
                H_[2,1] = H_[1,2]
                H_[1,3] = POWER*LOGRN^2
                H_[3,1] = H_[1,3]
                H_[2,2] = EV_[1]*POWER*LOGRN^2
                H_[2,3] = EV_[1]*POWER*LOGRN^3
                H_[3,2] = H_[2,3]
                H_[3,3] = EV_[1]*POWER*LOGRN^4
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

