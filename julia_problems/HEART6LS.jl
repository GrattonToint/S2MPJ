function HEART6LS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HEART6LS
#    *********
# 
#    Dipole model of the heart (6 x 6 version).
# 
#    This is a least-squares variant of problem HEART6.
# 
#    Source:
#    J. E. Dennis, Jr., D. M. Gay, P. A. Vu,
#    "A New Nonlinear Equations Test Problem".
#    Tech. Rep. 83-16, Dept. of Math. Sci., Rice Univ., Houston, TX
#    June 1983, revised May 1985.
# 
#    SIF input: A.R. Conn, May 1993, modified by Ph. Toint, March 1994.
#               correction by S. Gratton & Ph. Toint, May 2024
# 
#    classification = "C-SUR2-MN-6-0"
# 
# 
#    some useful parameters.
# 
# RE sum_Mx              0.485
# RE sum_My              -0.0019
# RE sum_A               -0.0581
# RE sum_B               0.015
# RE sum_C               0.105
# RE sum_D               0.0406
# RE sum_E               0.167
# RE sum_F               -0.399
# 
# RE sum_Mx              -0.69
# RE sum_My              -0.044
# RE sum_A               -1.57
# RE sum_B               -1.31
# RE sum_C               -2.65
# RE sum_D               2.0
# RE sum_E               -12.6
# RE sum_F               9.48
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HEART6LS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["sumuMx"] = -0.816
        v_["sumuMy"] = -0.017
        v_["sumuA"] = -1.826
        v_["sumuB"] = -0.754
        v_["sumuC"] = -4.839
        v_["sumuD"] = -3.259
        v_["sumuE"] = -14.023
        v_["sumuF"] = 15.467
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("a",ix_)
        arrset(pb.xnames,iv,"a")
        iv,ix_,_ = s2mpj_ii("c",ix_)
        arrset(pb.xnames,iv,"c")
        iv,ix_,_ = s2mpj_ii("t",ix_)
        arrset(pb.xnames,iv,"t")
        iv,ix_,_ = s2mpj_ii("u",ix_)
        arrset(pb.xnames,iv,"u")
        iv,ix_,_ = s2mpj_ii("v",ix_)
        arrset(pb.xnames,iv,"v")
        iv,ix_,_ = s2mpj_ii("w",ix_)
        arrset(pb.xnames,iv,"w")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("G1",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("G2",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("G3",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("G4",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("G5",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("G6",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["G1"]] = Float64(v_["sumuA"])
        pbm.gconst[ig_["G2"]] = Float64(v_["sumuB"])
        pbm.gconst[ig_["G3"]] = Float64(v_["sumuC"])
        pbm.gconst[ig_["G4"]] = Float64(v_["sumuD"])
        pbm.gconst[ig_["G5"]] = Float64(v_["sumuE"])
        pbm.gconst[ig_["G6"]] = Float64(v_["sumuF"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.x0[ix_["a"]] = Float64(0.0)
        pb.x0[ix_["c"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "en3PROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eVPV", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"ALPHA")
        it,iet_,_ = s2mpj_ii( "eADFSQ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "ePDFSQ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        loaset(elftp,it,1,"ALPHA")
        it,iet_,_ = s2mpj_ii( "eP3PRD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        loaset(elftp,it,1,"ALPHA")
        it,iet_,_ = s2mpj_ii( "en3DPRD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eD3PRD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        loaset(elftp,it,1,"ALPHA")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eVPV")
        arrset(ielftype,ie,iet_["eVPV"])
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eVPV")
        arrset(ielftype,ie,iet_["eVPV"])
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eVPV")
        arrset(ielftype,ie,iet_["eVPV"])
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
        ename = "E7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PROD")
        arrset(ielftype,ie,iet_["en2PROD"])
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eVPV")
        arrset(ielftype,ie,iet_["eVPV"])
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eADFSQ")
        arrset(ielftype,ie,iet_["eADFSQ"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PROD")
        arrset(ielftype,ie,iet_["en3PROD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePDFSQ")
        arrset(ielftype,ie,iet_["ePDFSQ"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
        ename = "E12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eP3PRD")
        arrset(ielftype,ie,iet_["eP3PRD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E13"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eADFSQ")
        arrset(ielftype,ie,iet_["eADFSQ"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E14"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PROD")
        arrset(ielftype,ie,iet_["en3PROD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E15"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePDFSQ")
        arrset(ielftype,ie,iet_["ePDFSQ"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E16"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eP3PRD")
        arrset(ielftype,ie,iet_["eP3PRD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
        ename = "E17"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3DPRD")
        arrset(ielftype,ie,iet_["en3DPRD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E18"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3DPRD")
        arrset(ielftype,ie,iet_["en3DPRD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E19"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eD3PRD")
        arrset(ielftype,ie,iet_["eD3PRD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
        ename = "E20"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eD3PRD")
        arrset(ielftype,ie,iet_["eD3PRD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E21"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3DPRD")
        arrset(ielftype,ie,iet_["en3DPRD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E22"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3DPRD")
        arrset(ielftype,ie,iet_["en3DPRD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "v"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "t"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E23"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eD3PRD")
        arrset(ielftype,ie,iet_["eD3PRD"])
        vname = "c"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMy"]))
        ename = "E24"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eD3PRD")
        arrset(ielftype,ie,iet_["eD3PRD"])
        vname = "a"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "w"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "u"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="ALPHA",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["sumuMx"]))
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
        ig = ig_["G1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["G2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E7"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E8"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["G3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E9"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E10"])
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E11"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E12"])
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        ig = ig_["G4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E13"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E14"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E15"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E16"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        ig = ig_["G5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E17"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E18"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E19"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E20"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["G6"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E21"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E22"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E23"])
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E24"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-MN-6-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en2PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = ZERO
                H_[2,2] = ZERO
                H_[1,2] = 1.0
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

    elseif action == "en3PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = ZERO
                H_[3,3] = ZERO
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
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

    elseif action == "eVPV"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        DIFF = pbm.elpar[iel_][1]-EV_[2]
        f_   = EV_[1]*DIFF
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DIFF
            g_[2] = -EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = ZERO
                H_[2,2] = ZERO
                H_[1,2] = -1.0
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

    elseif action == "eADFSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        TWO = 2.0
        DFSQ = EV_[2]^2-EV_[3]^2
        TWOX = TWO*EV_[1]
        f_   = EV_[1]*DFSQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DFSQ
            g_[2] = TWOX*EV_[2]
            g_[3] = -TWOX*EV_[3]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = TWOX
                H_[3,3] = -TWOX
                H_[1,2] = TWO*EV_[2]
                H_[2,1] = H_[1,2]
                H_[1,3] = -TWO*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = ZERO
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

    elseif action == "ePDFSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        TWO = 2.0
        DIFF = pbm.elpar[iel_][1]-EV_[1]
        DFSQ = EV_[2]^2-EV_[3]^2
        TWOX = TWO*EV_[1]
        TWOD = TWO*DIFF
        f_   = DIFF*DFSQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -DFSQ
            g_[2] = TWOD*EV_[2]
            g_[3] = -TWOD*EV_[3]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = TWOD
                H_[3,3] = -TWOD
                H_[1,2] = -TWO*EV_[2]
                H_[2,1] = H_[1,2]
                H_[1,3] = TWO*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = ZERO
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

    elseif action == "eP3PRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        TWO = 2.0
        DIFF = pbm.elpar[iel_][1]-EV_[1]
        f_   = DIFF*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EV_[2]*EV_[3]
            g_[2] = DIFF*EV_[3]
            g_[3] = DIFF*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = ZERO
                H_[3,3] = ZERO
                H_[1,2] = -EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = -EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = DIFF
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

    elseif action == "en3DPRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        TWO = 2.0
        THREE = 3.0
        SIX = 6.0
        DIFF = EV_[2]^2-THREE*EV_[3]^2
        f_   = EV_[1]*EV_[2]*DIFF
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*DIFF
            g_[2] = EV_[1]*DIFF+TWO*EV_[1]*EV_[2]^2
            g_[3] = -SIX*EV_[1]*EV_[2]*EV_[3]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = SIX*EV_[1]*EV_[2]
                H_[3,3] = -SIX*EV_[1]*EV_[2]
                H_[1,2] = DIFF+TWO*EV_[2]^2
                H_[2,1] = H_[1,2]
                H_[1,3] = -SIX*EV_[2]*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = -SIX*EV_[1]*EV_[3]
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

    elseif action == "eD3PRD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        ZERO = 0.0
        TWO = 2.0
        THREE = 3.0
        SIX = 6.0
        DFSQ = EV_[2]^2-THREE*EV_[3]^2
        DIFF = pbm.elpar[iel_][1]-EV_[1]
        f_   = DIFF*EV_[2]*DFSQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EV_[2]*DFSQ
            g_[2] = DIFF*(DFSQ+TWO*EV_[2]^2)
            g_[3] = -SIX*EV_[2]*EV_[3]*DIFF
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = ZERO
                H_[2,2] = SIX*EV_[2]*DIFF
                H_[3,3] = -SIX*EV_[2]*DIFF
                H_[1,2] = -DFSQ-TWO*EV_[2]^2
                H_[2,1] = H_[1,2]
                H_[1,3] = SIX*EV_[2]*EV_[3]
                H_[3,1] = H_[1,3]
                H_[2,3] = -SIX*DIFF*EV_[3]
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

