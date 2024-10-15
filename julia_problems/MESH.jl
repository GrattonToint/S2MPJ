function MESH(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    The goodness of a finite element grid is characterized by the
#    smallest angle of all triangles. Given a triangulation of a domain 
#    in R**2. Find a topological equivalent triangulation so, that the 
#    smallest angel becomes as large as possible. Topological equivalent 
#    means shifting the edges of the grid only in such a way that 
#    neighbouring triangles remain neighbours. 
# 
#    Source: Prof. Dr. Michael Kraetzschmar, Institut fuer Angewandte 
#            Mathematik der Fachhochschule Flensburg, Kanzleistrasse 91-93, 
#            D-24943 FLENSBURG, GERMANY
# 
#    SIF input: Prof. Dr. Michael Kraetzschmar
# 
#    classification = "C-OOR2-AY-41-48"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MESH"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["omega1"] = 1.0e3
        v_["omega2"] = -1.0e3
        v_["omega3"] = -1.0e5
        v_["s"] = 0.700000
        v_["pi"] = acos(-1.0)
        v_["sqrt3/2"] = sqrt(0.75)
        v_["h"] = v_["sqrt3/2"]*v_["s"]
        v_["drei"] = 3.0
        v_["pi/3"] = v_["pi"]/v_["drei"]
        v_["1"] = 1
        v_["np"] = 5
        v_["nk"] = 8
        v_["nd"] = 4
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for i = Int64(v_["1"]):Int64(v_["np"])
            iv,ix_,_ = s2mpj_ii("x"*string(i),ix_)
            arrset(pb.xnames,iv,"x"*string(i))
            iv,ix_,_ = s2mpj_ii("y"*string(i),ix_)
            arrset(pb.xnames,iv,"y"*string(i))
        end
        for i = Int64(v_["1"]):Int64(v_["nk"])
            iv,ix_,_ = s2mpj_ii("l"*string(i),ix_)
            arrset(pb.xnames,iv,"l"*string(i))
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            iv,ix_,_ = s2mpj_ii("alpha"*string(i),ix_)
            arrset(pb.xnames,iv,"alpha"*string(i))
            iv,ix_,_ = s2mpj_ii("beta"*string(i),ix_)
            arrset(pb.xnames,iv,"beta"*string(i))
            iv,ix_,_ = s2mpj_ii("gamma"*string(i),ix_)
            arrset(pb.xnames,iv,"gamma"*string(i))
            iv,ix_,_ = s2mpj_ii("delta"*string(i),ix_)
            arrset(pb.xnames,iv,"delta"*string(i))
            iv,ix_,_ = s2mpj_ii("f"*string(i),ix_)
            arrset(pb.xnames,iv,"f"*string(i))
        end
        iv,ix_,_ = s2mpj_ii("deltamin",ix_)
        arrset(pb.xnames,iv,"deltamin")
        iv,ix_,_ = s2mpj_ii("fmin",ix_)
        arrset(pb.xnames,iv,"fmin")
        iv,ix_,_ = s2mpj_ii("fmax",ix_)
        arrset(pb.xnames,iv,"fmax")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("obj1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["deltamin"]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(v_["omega1"]))
        ig,ig_,_ = s2mpj_ii("obj2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["fmax"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["fmin"]
        pbm.A[ig,iv] += Float64(-1.0)
        arrset(pbm.gscale,ig,Float64(v_["omega2"]))
        ig,ig_,_ = s2mpj_ii("obj3",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["omega3"]))
        for i = Int64(v_["1"]):Int64(v_["nk"])
            ig,ig_,_ = s2mpj_ii("seit"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"seit"*string(i))
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig,ig_,_ = s2mpj_ii("skal"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"skal"*string(i))
            ig,ig_,_ = s2mpj_ii("skbe"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"skbe"*string(i))
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig,ig_,_ = s2mpj_ii("doppf"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"doppf"*string(i))
            iv = ix_["f"*string(i)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig,ig_,_ = s2mpj_ii("wisum"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"wisum"*string(i))
            iv = ix_["alpha"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["beta"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["gamma"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig,ig_,_ = s2mpj_ii("alphd"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"alphd"*string(i))
            iv = ix_["alpha"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["delta"*string(i)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("betad"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"betad"*string(i))
            iv = ix_["beta"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["delta"*string(i)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("gammd"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"gammd"*string(i))
            iv = ix_["gamma"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["delta"*string(i)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("deltd"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"deltd"*string(i))
            iv = ix_["delta"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["deltamin"]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig,ig_,_ = s2mpj_ii("fmind"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"fmind"*string(i))
            iv = ix_["f"*string(i)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["fmin"]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("fmaxd"*string(i),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"fmaxd"*string(i))
            iv = ix_["fmax"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["f"*string(i)]
            pbm.A[ig,iv] += Float64(-1.0)
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
        for i = Int64(v_["1"]):Int64(v_["nd"])
            pbm.gconst[ig_["wisum"*string(i)]] = Float64(v_["pi"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["deltamin"]] = v_["pi"]
        pb.xlower[ix_["x1"]] = 0.000000
        pb.xupper[ix_["x1"]] = 0.000000
        pb.xlower[ix_["y1"]] = 0.000000
        pb.xupper[ix_["y1"]] = 0.000000
        pb.xlower[ix_["x2"]] = 0.000000
        pb.xupper[ix_["x2"]] = 0.000000
        pb.xlower[ix_["y2"]] = 1.000000
        pb.xupper[ix_["y2"]] = 1.000000
        pb.xlower[ix_["x3"]] = 1.000000
        pb.xupper[ix_["x3"]] = 1.000000
        pb.xlower[ix_["y3"]] = 1.000000
        pb.xupper[ix_["y3"]] = 1.000000
        pb.xlower[ix_["x4"]] = 1.000000
        pb.xupper[ix_["x4"]] = 1.000000
        pb.xlower[ix_["y4"]] = 0.000000
        pb.xupper[ix_["y4"]] = 0.000000
        pb.xlower[ix_["x5"]] = -Inf
        pb.xupper[ix_["x5"]] = +Inf
        pb.xlower[ix_["y5"]] = -Inf
        pb.xupper[ix_["y5"]] = +Inf
        for i = Int64(v_["1"]):Int64(v_["nd"])
            pb.xupper[ix_["alpha"*string(i)]] = v_["pi"]
            pb.xupper[ix_["beta"*string(i)]] = v_["pi"]
            pb.xupper[ix_["gamma"*string(i)]] = v_["pi"]
            pb.xupper[ix_["delta"*string(i)]] = v_["pi"]
        end
        pb.xupper[ix_["deltamin"]] = v_["pi"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["x1"]] = Float64(0.000000)
        pb.x0[ix_["y1"]] = Float64(0.000000)
        pb.x0[ix_["x2"]] = Float64(0.000000)
        pb.x0[ix_["y2"]] = Float64(1.000000)
        pb.x0[ix_["x3"]] = Float64(1.000000)
        pb.x0[ix_["y3"]] = Float64(1.000000)
        pb.x0[ix_["x4"]] = Float64(1.000000)
        pb.x0[ix_["y4"]] = Float64(0.000000)
        pb.x0[ix_["x5"]] = Float64(0.350000)
        pb.x0[ix_["y5"]] = Float64(0.606218)
        pb.x0[ix_["l1"]] = Float64(0.700000)
        pb.x0[ix_["l2"]] = Float64(0.526844)
        pb.x0[ix_["l3"]] = Float64(1.000000)
        pb.x0[ix_["l4"]] = Float64(0.759977)
        pb.x0[ix_["l5"]] = Float64(1.000000)
        pb.x0[ix_["l6"]] = Float64(0.888819)
        pb.x0[ix_["l7"]] = Float64(1.000000)
        pb.x0[ix_["l8"]] = Float64(1.000000)
        pb.x0[ix_["alpha1"]] = Float64(0.523599)
        pb.x0[ix_["beta1"]] = Float64(1.891392)
        pb.x0[ix_["gamma1"]] = Float64(0.726602)
        pb.x0[ix_["delta1"]] = Float64(0.523599)
        pb.x0[ix_["f1"]] = Float64(0.350000)
        pb.x0[ix_["alpha2"]] = Float64(0.844195)
        pb.x0[ix_["beta2"]] = Float64(1.752711)
        pb.x0[ix_["gamma2"]] = Float64(0.544687)
        pb.x0[ix_["delta2"]] = Float64(0.544687)
        pb.x0[ix_["f2"]] = Float64(0.393782)
        pb.x0[ix_["alpha3"]] = Float64(1.026109)
        pb.x0[ix_["beta3"]] = Float64(1.295247)
        pb.x0[ix_["gamma3"]] = Float64(0.820236)
        pb.x0[ix_["delta3"]] = Float64(0.820236)
        pb.x0[ix_["f3"]] = Float64(0.650000)
        pb.x0[ix_["alpha4"]] = Float64(0.750560)
        pb.x0[ix_["beta4"]] = Float64(1.343835)
        pb.x0[ix_["gamma4"]] = Float64(1.047198)
        pb.x0[ix_["delta4"]] = Float64(0.750560)
        pb.x0[ix_["f4"]] = Float64(0.606218)
        pb.x0[ix_["deltamin"]] = Float64(0.523599)
        pb.x0[ix_["fmin"]] = Float64(0.350000)
        pb.x0[ix_["fmax"]] = Float64(0.650000)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ediffsq", iet_)
        loaset(elftv,it,1,"winkel")
        loaset(elftv,it,2,"minwi")
        it,iet_,_ = s2mpj_ii( "elaenge", iet_)
        loaset(elftv,it,1,"lp1x")
        loaset(elftv,it,2,"lp1y")
        loaset(elftv,it,3,"lp2x")
        loaset(elftv,it,4,"lp2y")
        loaset(elftv,it,5,"llaenge")
        it,iet_,_ = s2mpj_ii( "evekprod", iet_)
        loaset(elftv,it,1,"vp1x")
        loaset(elftv,it,2,"vp1y")
        loaset(elftv,it,3,"vp2x")
        loaset(elftv,it,4,"vp2y")
        loaset(elftv,it,5,"vp3x")
        loaset(elftv,it,6,"vp3y")
        it,iet_,_ = s2mpj_ii( "esklprod", iet_)
        loaset(elftv,it,1,"sp1x")
        loaset(elftv,it,2,"sp1y")
        loaset(elftv,it,3,"sp2x")
        loaset(elftv,it,4,"sp2y")
        loaset(elftv,it,5,"sp3x")
        loaset(elftv,it,6,"sp3y")
        it,iet_,_ = s2mpj_ii( "ecosprod", iet_)
        loaset(elftv,it,1,"cl1")
        loaset(elftv,it,2,"cl2")
        loaset(elftv,it,3,"cw")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "aldsq"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ediffsq")
            arrset(ielftype,ie,iet_["ediffsq"])
            vname = "alpha"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="winkel",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "delta"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="minwi",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "bedsq"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ediffsq")
            arrset(ielftype,ie,iet_["ediffsq"])
            vname = "beta"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="winkel",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "delta"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="minwi",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "gadsq"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ediffsq")
            arrset(ielftype,ie,iet_["ediffsq"])
            vname = "gamma"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="winkel",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "delta"*string(i)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="minwi",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for i = Int64(v_["1"]):Int64(v_["nk"])
            ename = "laeng"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"elaenge")
            arrset(ielftype,ie,iet_["elaenge"])
        end
        ename = "laeng1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "laeng8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="lp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="llaenge",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "sal"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"esklprod")
            arrset(ielftype,ie,iet_["esklprod"])
        end
        ename = "sal1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sal2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sal3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sal4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "cal"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ecosprod")
            arrset(ielftype,ie,iet_["ecosprod"])
        end
        ename = "cal1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "alpha1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cal2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "alpha2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cal3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "alpha3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cal4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "alpha4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "sbe"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"esklprod")
            arrset(ielftype,ie,iet_["esklprod"])
        end
        ename = "sbe1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sbe2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sbe3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "sbe4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="sp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "cbe"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ecosprod")
            arrset(ielftype,ie,iet_["ecosprod"])
        end
        ename = "cbe1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "beta1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cbe2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "beta2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cbe3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "beta3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "cbe4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "l6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "l1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cl2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "beta4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="cw",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ename = "flae"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"evekprod")
            arrset(ielftype,ie,iet_["evekprod"])
        end
        ename = "flae1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "flae2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "flae3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "flae4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "x4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp1y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp2y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "x1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3x",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "y1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="vp3y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gsquare",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["obj1"]
        arrset(pbm.grftype,ig,"gsquare")
        ig = ig_["obj2"]
        arrset(pbm.grftype,ig,"gsquare")
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig = ig_["obj3"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["aldsq"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["bedsq"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["gadsq"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        for i = Int64(v_["1"]):Int64(v_["nk"])
            ig = ig_["seit"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["laeng"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig = ig_["skal"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["sal"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["cal"*string(i)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["skbe"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["sbe"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["cbe"*string(i)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        for i = Int64(v_["1"]):Int64(v_["nd"])
            ig = ig_["doppf"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["flae"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              5.9213448D-4
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AY-41-48"
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

    elseif action == "ediffsq"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "elaenge"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,5)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        U_[3,5] = U_[3,5]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[3]*IV_[3]-IV_[1]*IV_[1]-IV_[2]*IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[3] = IV_[3]+IV_[3]
            g_[1] = -(IV_[1]+IV_[1])
            g_[2] = -(IV_[2]+IV_[2])
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[3,3] = 2.0
                H_[3,1] = 0.0
                H_[1,3] = H_[3,1]
                H_[3,2] = 0.0
                H_[2,3] = H_[3,2]
                H_[1,2] = 0.0
                H_[2,1] = H_[1,2]
                H_[1,1] = -2.0
                H_[2,2] = -2.0
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "evekprod"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,4,6)
        IV_ =  zeros(Float64,4)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        U_[3,5] = U_[3,5]+1
        U_[3,3] = U_[3,3]-1
        U_[4,6] = U_[4,6]+1
        U_[4,4] = U_[4,4]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        IV_[4] = dot(U_[4,:],EV_)
        f_   = IV_[3]*IV_[2]-IV_[1]*IV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -IV_[4]
            g_[2] = IV_[3]
            g_[3] = IV_[2]
            g_[4] = -IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[3,2] = 1.0
                H_[2,3] = H_[3,2]
                H_[1,4] = -1.0
                H_[4,1] = H_[1,4]
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "esklprod"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,4,6)
        IV_ =  zeros(Float64,4)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        U_[3,5] = U_[3,5]+1
        U_[3,3] = U_[3,3]-1
        U_[4,6] = U_[4,6]+1
        U_[4,4] = U_[4,4]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        IV_[4] = dot(U_[4,:],EV_)
        f_   = IV_[1]*IV_[3]+IV_[2]*IV_[4]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[3]
            g_[2] = IV_[4]
            g_[3] = IV_[1]
            g_[4] = IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,3] = 1.0
                H_[3,1] = H_[1,3]
                H_[2,4] = 1.0
                H_[4,2] = H_[2,4]
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ecosprod"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        cosa = cos(EV_[3])
        sina = sin(EV_[3])
        prod2 = EV_[1]*EV_[2]
        prod3 = prod2*cosa
        f_   = prod3
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*cosa
            g_[2] = EV_[1]*cosa
            g_[3] = -sina*prod2
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = cosa
                H_[2,1] = H_[1,2]
                H_[1,3] = -sina*EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = -sina*EV_[1]
                H_[3,2] = H_[2,3]
                H_[3,3] = -prod3
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

    elseif action == "gsquare"

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

