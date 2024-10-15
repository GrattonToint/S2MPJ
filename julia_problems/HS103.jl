function HS103(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: problem 103 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: N. Gould, December 1989.
# 
#    classification = "C-OOR2-AN-7-5"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS103"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["M"] = 5
        v_["N"] = 7
        v_["A101"] = 0.5
        v_["A102"] = 0.125
        v_["A103"] = 0.5
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("CONSTR"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"CONSTR"*string(I))
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
        pbm.gconst[ig_["CONSTR1"]] = Float64(1.0)
        pbm.gconst[ig_["CONSTR2"]] = Float64(1.0)
        pbm.gconst[ig_["CONSTR3"]] = Float64(1.0)
        pbm.gconst[ig_["CONSTR4"]] = Float64(1.0)
        pbm.gconst[ig_["CONSTR5"]] = Float64(3000.0)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.1,pb.n)
        pb.xupper = fill(10.0,pb.n)
        pb.xlower[ix_["X7"]] = 0.01
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(6.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en3PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        it,iet_,_ = s2mpj_ii( "en4PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        it,iet_,_ = s2mpj_ii( "en5PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        loaset(elftp,it,5,"P5")
        it,iet_,_ = s2mpj_ii( "en6PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P6")
        loaset(elftp,it,3,"P2")
        loaset(elftp,it,4,"P3")
        loaset(elftp,it,5,"P4")
        loaset(elftp,it,6,"P5")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1C1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en4PR")
        arrset(ielftype,ie,iet_["en4PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E2C1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        ename = "E3C1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.5))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.66666666))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.25))
        ename = "E1C2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.5))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E2C2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en4PR")
        arrset(ielftype,ie,iet_["en4PR"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        ename = "E3C2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.3333333333))
        ename = "E1C3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.5))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.3333333333))
        ename = "E2C3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.5))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.5))
        ename = "E3C3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en4PR")
        arrset(ielftype,ie,iet_["en4PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E4C3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E1C4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.3333333333))
        ename = "E2C4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en6PR")
        arrset(ielftype,ie,iet_["en6PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.3333333333))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.666666666))
        posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.25))
        ename = "E3C4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-3.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.75))
        ename = "E4C4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PR")
        arrset(ielftype,ie,iet_["en3PR"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        ename = "E1C5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-3.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["A101"]))
        ename = "E2C5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en6PR")
        arrset(ielftype,ie,iet_["en6PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.5))
        ename = "E3C5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en5PR")
        arrset(ielftype,ie,iet_["en5PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "E4C5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en6PR")
        arrset(ielftype,ie,iet_["en6PR"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),Float64(10.0),Float64(6.0)))
        posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.5))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
        posep = findfirst(x->x=="P6",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(10.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C5"])
        loaset(pbm.grelw,ig,posel,Float64(15.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(20.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4C5"])
        loaset(pbm.grelw,ig,posel,Float64(25.0))
        ig = ig_["CONSTR1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C1"])
        loaset(pbm.grelw,ig,posel,Float64(0.7))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.2))
        ig = ig_["CONSTR2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.3))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C2"])
        loaset(pbm.grelw,ig,posel,Float64(0.8))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(3.1))
        ig = ig_["CONSTR3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C3"])
        loaset(pbm.grelw,ig,posel,Float64(0.1))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4C3"])
        loaset(pbm.grelw,ig,posel,Float64(0.65))
        ig = ig_["CONSTR4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.2))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C4"])
        loaset(pbm.grelw,ig,posel,Float64(0.3))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.4))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4C4"])
        loaset(pbm.grelw,ig,posel,Float64(0.5))
        ig = ig_["CONSTR5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(10.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2C5"])
        loaset(pbm.grelw,ig,posel,Float64(15.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3C5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(20.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4C5"])
        loaset(pbm.grelw,ig,posel,Float64(25.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               1809.76476
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[1:pb.nle] = grange[legrps]
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AN-7-5"
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

    elseif action == "en3PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVALUE  = (
              (EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^pbm.elpar[iel_][2])*(EV_[3]^pbm.elpar[iel_][3]))
        f_   = FVALUE
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[3])*((pbm.elpar[iel_][3]-1.0)/EV_[3]))
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,2] = H_[2,3]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en4PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVALUE  = (
              (EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^pbm.elpar[iel_][2])*(EV_[3]^pbm.elpar[iel_][3])*(EV_[4]^pbm.elpar[iel_][4]))
        f_   = FVALUE
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])
            g_[4] = FVALUE*(pbm.elpar[iel_][4]/EV_[4])
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[3])*((pbm.elpar[iel_][3]-1.0)/EV_[3]))
                H_[4,4]  = (
                      FVALUE*(pbm.elpar[iel_][4]/EV_[4])*((pbm.elpar[iel_][4]-1.0)/EV_[4]))
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,1] = H_[1,4]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,2] = H_[2,4]
                H_[3,4] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,3] = H_[3,4]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en5PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVALUE  = (
              (EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^pbm.elpar[iel_][2])*(EV_[3]^pbm.elpar[iel_][3])*(EV_[4]^pbm.elpar[iel_][4])*(EV_[5]^pbm.elpar[iel_][5]))
        f_   = FVALUE
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])
            g_[4] = FVALUE*(pbm.elpar[iel_][4]/EV_[4])
            g_[5] = FVALUE*(pbm.elpar[iel_][5]/EV_[5])
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[2])*((pbm.elpar[iel_][2]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[3])*((pbm.elpar[iel_][3]-1.0)/EV_[3]))
                H_[4,4]  = (
                      FVALUE*(pbm.elpar[iel_][4]/EV_[4])*((pbm.elpar[iel_][4]-1.0)/EV_[4]))
                H_[5,5]  = (
                      FVALUE*(pbm.elpar[iel_][5]/EV_[5])*((pbm.elpar[iel_][5]-1.0)/EV_[5]))
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,1] = H_[1,4]
                H_[1,5] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][5]/EV_[5])
                H_[5,1] = H_[1,5]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][3]/EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,2] = H_[2,4]
                H_[2,5] = FVALUE*(pbm.elpar[iel_][2]/EV_[2])*(pbm.elpar[iel_][5]/EV_[5])
                H_[5,2] = H_[2,5]
                H_[3,4] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])*(pbm.elpar[iel_][4]/EV_[4])
                H_[4,3] = H_[3,4]
                H_[3,5] = FVALUE*(pbm.elpar[iel_][3]/EV_[3])*(pbm.elpar[iel_][5]/EV_[5])
                H_[5,3] = H_[3,5]
                H_[4,5] = FVALUE*(pbm.elpar[iel_][4]/EV_[4])*(pbm.elpar[iel_][5]/EV_[5])
                H_[5,4] = H_[4,5]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en6PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        FVALUE  = (
              (EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^pbm.elpar[iel_][3])*(EV_[3]^pbm.elpar[iel_][4])*(EV_[4]^pbm.elpar[iel_][5])*(EV_[5]^pbm.elpar[iel_][6])*(EV_[6]^pbm.elpar[iel_][2]))
        f_   = FVALUE
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])
            g_[2] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])
            g_[3] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])
            g_[4] = FVALUE*(pbm.elpar[iel_][5]/EV_[4])
            g_[5] = FVALUE*(pbm.elpar[iel_][6]/EV_[5])
            g_[6] = FVALUE*(pbm.elpar[iel_][2]/EV_[6])
            if nargout>2
                H_ = zeros(Float64,6,6)
                H_[1,1]  = (
                      FVALUE*(pbm.elpar[iel_][1]/EV_[1])*((pbm.elpar[iel_][1]-1.0)/EV_[1]))
                H_[2,2]  = (
                      FVALUE*(pbm.elpar[iel_][3]/EV_[2])*((pbm.elpar[iel_][3]-1.0)/EV_[2]))
                H_[3,3]  = (
                      FVALUE*(pbm.elpar[iel_][4]/EV_[3])*((pbm.elpar[iel_][4]-1.0)/EV_[3]))
                H_[4,4]  = (
                      FVALUE*(pbm.elpar[iel_][5]/EV_[4])*((pbm.elpar[iel_][5]-1.0)/EV_[4]))
                H_[5,5]  = (
                      FVALUE*(pbm.elpar[iel_][6]/EV_[5])*((pbm.elpar[iel_][6]-1.0)/EV_[5]))
                H_[6,6]  = (
                      FVALUE*(pbm.elpar[iel_][2]/EV_[6])*((pbm.elpar[iel_][2]-1.0)/EV_[6]))
                H_[1,2] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][3]/EV_[2])
                H_[2,1] = H_[1,2]
                H_[1,3] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][4]/EV_[3])
                H_[3,1] = H_[1,3]
                H_[1,4] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,1] = H_[1,4]
                H_[1,5] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][6]/EV_[5])
                H_[5,1] = H_[1,5]
                H_[1,6] = FVALUE*(pbm.elpar[iel_][1]/EV_[1])*(pbm.elpar[iel_][2]/EV_[6])
                H_[6,1] = H_[1,6]
                H_[2,3] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][4]/EV_[3])
                H_[3,2] = H_[2,3]
                H_[2,4] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,2] = H_[2,4]
                H_[2,5] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][6]/EV_[5])
                H_[5,2] = H_[2,5]
                H_[2,6] = FVALUE*(pbm.elpar[iel_][3]/EV_[2])*(pbm.elpar[iel_][2]/EV_[6])
                H_[6,2] = H_[2,6]
                H_[3,4] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])*(pbm.elpar[iel_][5]/EV_[4])
                H_[4,3] = H_[3,4]
                H_[3,5] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])*(pbm.elpar[iel_][6]/EV_[5])
                H_[5,3] = H_[3,5]
                H_[3,6] = FVALUE*(pbm.elpar[iel_][4]/EV_[3])*(pbm.elpar[iel_][2]/EV_[6])
                H_[6,3] = H_[3,6]
                H_[4,5] = FVALUE*(pbm.elpar[iel_][5]/EV_[4])*(pbm.elpar[iel_][6]/EV_[5])
                H_[5,4] = H_[4,5]
                H_[4,6] = FVALUE*(pbm.elpar[iel_][5]/EV_[4])*(pbm.elpar[iel_][2]/EV_[6])
                H_[6,4] = H_[4,6]
                H_[5,6] = FVALUE*(pbm.elpar[iel_][6]/EV_[5])*(pbm.elpar[iel_][2]/EV_[6])
                H_[6,5] = H_[5,6]
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

