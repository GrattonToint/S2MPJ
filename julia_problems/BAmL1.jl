function BAmL1(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BAmL1
#    *********
# 
#    Bundle Adjustment problem from reconstructive geometry in which
#    a collection of photographs is used to determine the position of
#    a set of observed points. Each observed point is seen via its
#    two-dimensional projections on a subset of the photographs. The
#    solution is found by solvng a large nonlinear least-squares problem.
#    This variant is given as an inconsistent set of nonlinear equations.
# 
#    Source: data from the Bundle Adjustment in the Large
#    project, http://grail.cs.washington.edu/projects/bal/
# 
#    Ladybug datasets (single image extracted)
# 
#    SIF input: Nick Gould, November 2016
# 
#    classification = "NOR2-MN-57-12"
# 
#    Number of images
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BAmL1"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "BAmL1"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["nuimages"] = 49
        v_["nupoints"] = 1
        v_["nuobservs"] = 6
        v_["1"] = 1
        v_["O1"] = 1.0
        v_["O2"] = 2.0
        v_["O3"] = 4.0
        v_["O4"] = 27.0
        v_["O5"] = 30.0
        v_["O6"] = 37.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["nupoints"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2x_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2x_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        for J = Int64(v_["1"]):Int64(v_["nuobservs"])
            v_["RI"] = v_["O"*string(J)]
            v_["I"] = trunc(Int,v_["RI"])
            iv,ix_,_ = s2x_ii("RX"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"RX"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("RY"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"RY"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("RZ"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"RZ"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("TX"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"TX"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("TY"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"TY"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("TZ"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"TZ"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("KA"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"KA"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("KB"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"KB"*string(Int64(v_["I"])))
            iv,ix_,_ = s2x_ii("F"*string(Int64(v_["I"])),ix_)
            arrset(pb.xnames,iv,"F"*string(Int64(v_["I"])))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["nuobservs"])
            ig,ig_,_ = s2x_ii("RX"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"RX"*string(I))
            ig,ig_,_ = s2x_ii("RY"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"RY"*string(I))
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["RX1"]] = Float64(-332.65)
        pbm.gconst[ig_["RY1"]] = Float64(262.09)
        pbm.gconst[ig_["RX2"]] = Float64(-199.76)
        pbm.gconst[ig_["RY2"]] = Float64(166.7)
        pbm.gconst[ig_["RX3"]] = Float64(-253.06)
        pbm.gconst[ig_["RY3"]] = Float64(202.27)
        pbm.gconst[ig_["RX4"]] = Float64(58.13)
        pbm.gconst[ig_["RY4"]] = Float64(271.89)
        pbm.gconst[ig_["RX5"]] = Float64(238.22)
        pbm.gconst[ig_["RY5"]] = Float64(237.37)
        pbm.gconst[ig_["RX6"]] = Float64(317.55)
        pbm.gconst[ig_["RY6"]] = Float64(221.15)
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(-.6120001572)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(-.6120001572)
        end
        if haskey(ix_,"Y1")
            pb.x0[ix_["Y1"]] = Float64(.57175904776)
        else
            pb.y0[findfirst(x->x==ig_["Y1"],pbm.congrps)] = Float64(.57175904776)
        end
        if haskey(ix_,"Z1")
            pb.x0[ix_["Z1"]] = Float64(-1.847081276)
        else
            pb.y0[findfirst(x->x==ig_["Z1"],pbm.congrps)] = Float64(-1.847081276)
        end
        if haskey(ix_,"RX1")
            pb.x0[ix_["RX1"]] = Float64(.01574151594)
        else
            pb.y0[findfirst(x->x==ig_["RX1"],pbm.congrps)] = Float64(.01574151594)
        end
        if haskey(ix_,"RY1")
            pb.x0[ix_["RY1"]] = Float64(-.0127909362)
        else
            pb.y0[findfirst(x->x==ig_["RY1"],pbm.congrps)] = Float64(-.0127909362)
        end
        if haskey(ix_,"RZ1")
            pb.x0[ix_["RZ1"]] = Float64(-.0044008498)
        else
            pb.y0[findfirst(x->x==ig_["RZ1"],pbm.congrps)] = Float64(-.0044008498)
        end
        if haskey(ix_,"TX1")
            pb.x0[ix_["TX1"]] = Float64(-.0340938396)
        else
            pb.y0[findfirst(x->x==ig_["TX1"],pbm.congrps)] = Float64(-.0340938396)
        end
        if haskey(ix_,"TY1")
            pb.x0[ix_["TY1"]] = Float64(-.107513871)
        else
            pb.y0[findfirst(x->x==ig_["TY1"],pbm.congrps)] = Float64(-.107513871)
        end
        if haskey(ix_,"TZ1")
            pb.x0[ix_["TZ1"]] = Float64(1.1202240291)
        else
            pb.y0[findfirst(x->x==ig_["TZ1"],pbm.congrps)] = Float64(1.1202240291)
        end
        if haskey(ix_,"KA1")
            pb.x0[ix_["KA1"]] = Float64(-3.177064E-7)
        else
            pb.y0[findfirst(x->x==ig_["KA1"],pbm.congrps)] = Float64(-3.177064E-7)
        end
        if haskey(ix_,"KB1")
            pb.x0[ix_["KB1"]] = Float64(5.882049E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB1"],pbm.congrps)] = Float64(5.882049E-13)
        end
        if haskey(ix_,"F1")
            pb.x0[ix_["F1"]] = Float64(399.75152639)
        else
            pb.y0[findfirst(x->x==ig_["F1"],pbm.congrps)] = Float64(399.75152639)
        end
        if haskey(ix_,"RX2")
            pb.x0[ix_["RX2"]] = Float64(.01597732412)
        else
            pb.y0[findfirst(x->x==ig_["RX2"],pbm.congrps)] = Float64(.01597732412)
        end
        if haskey(ix_,"RY2")
            pb.x0[ix_["RY2"]] = Float64(-.0252244646)
        else
            pb.y0[findfirst(x->x==ig_["RY2"],pbm.congrps)] = Float64(-.0252244646)
        end
        if haskey(ix_,"RZ2")
            pb.x0[ix_["RZ2"]] = Float64(-.0094001416)
        else
            pb.y0[findfirst(x->x==ig_["RZ2"],pbm.congrps)] = Float64(-.0094001416)
        end
        if haskey(ix_,"TX2")
            pb.x0[ix_["TX2"]] = Float64(-.0085667661)
        else
            pb.y0[findfirst(x->x==ig_["TX2"],pbm.congrps)] = Float64(-.0085667661)
        end
        if haskey(ix_,"TY2")
            pb.x0[ix_["TY2"]] = Float64(-.1218804907)
        else
            pb.y0[findfirst(x->x==ig_["TY2"],pbm.congrps)] = Float64(-.1218804907)
        end
        if haskey(ix_,"TZ2")
            pb.x0[ix_["TZ2"]] = Float64(.7190133075)
        else
            pb.y0[findfirst(x->x==ig_["TZ2"],pbm.congrps)] = Float64(.7190133075)
        end
        if haskey(ix_,"KA2")
            pb.x0[ix_["KA2"]] = Float64(-3.780477E-7)
        else
            pb.y0[findfirst(x->x==ig_["KA2"],pbm.congrps)] = Float64(-3.780477E-7)
        end
        if haskey(ix_,"KB2")
            pb.x0[ix_["KB2"]] = Float64(9.307431E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB2"],pbm.congrps)] = Float64(9.307431E-13)
        end
        if haskey(ix_,"F2")
            pb.x0[ix_["F2"]] = Float64(402.01753386)
        else
            pb.y0[findfirst(x->x==ig_["F2"],pbm.congrps)] = Float64(402.01753386)
        end
        if haskey(ix_,"RX4")
            pb.x0[ix_["RX4"]] = Float64(.01484625118)
        else
            pb.y0[findfirst(x->x==ig_["RX4"],pbm.congrps)] = Float64(.01484625118)
        end
        if haskey(ix_,"RY4")
            pb.x0[ix_["RY4"]] = Float64(-.0210628994)
        else
            pb.y0[findfirst(x->x==ig_["RY4"],pbm.congrps)] = Float64(-.0210628994)
        end
        if haskey(ix_,"RZ4")
            pb.x0[ix_["RZ4"]] = Float64(-.001166948)
        else
            pb.y0[findfirst(x->x==ig_["RZ4"],pbm.congrps)] = Float64(-.001166948)
        end
        if haskey(ix_,"TX4")
            pb.x0[ix_["TX4"]] = Float64(-.0249509707)
        else
            pb.y0[findfirst(x->x==ig_["TX4"],pbm.congrps)] = Float64(-.0249509707)
        end
        if haskey(ix_,"TY4")
            pb.x0[ix_["TY4"]] = Float64(-.1139847055)
        else
            pb.y0[findfirst(x->x==ig_["TY4"],pbm.congrps)] = Float64(-.1139847055)
        end
        if haskey(ix_,"TZ4")
            pb.x0[ix_["TZ4"]] = Float64(.92166020737)
        else
            pb.y0[findfirst(x->x==ig_["TZ4"],pbm.congrps)] = Float64(.92166020737)
        end
        if haskey(ix_,"KA4")
            pb.x0[ix_["KA4"]] = Float64(-3.295265E-7)
        else
            pb.y0[findfirst(x->x==ig_["KA4"],pbm.congrps)] = Float64(-3.295265E-7)
        end
        if haskey(ix_,"KB4")
            pb.x0[ix_["KB4"]] = Float64(6.732885E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB4"],pbm.congrps)] = Float64(6.732885E-13)
        end
        if haskey(ix_,"F4")
            pb.x0[ix_["F4"]] = Float64(400.40175368)
        else
            pb.y0[findfirst(x->x==ig_["F4"],pbm.congrps)] = Float64(400.40175368)
        end
        if haskey(ix_,"RX27")
            pb.x0[ix_["RX27"]] = Float64(.01991666998)
        else
            pb.y0[findfirst(x->x==ig_["RX27"],pbm.congrps)] = Float64(.01991666998)
        end
        if haskey(ix_,"RY27")
            pb.x0[ix_["RY27"]] = Float64(-1.22433082)
        else
            pb.y0[findfirst(x->x==ig_["RY27"],pbm.congrps)] = Float64(-1.22433082)
        end
        if haskey(ix_,"RZ27")
            pb.x0[ix_["RZ27"]] = Float64(.0119988756)
        else
            pb.y0[findfirst(x->x==ig_["RZ27"],pbm.congrps)] = Float64(.0119988756)
        end
        if haskey(ix_,"TX27")
            pb.x0[ix_["TX27"]] = Float64(-1.411897512)
        else
            pb.y0[findfirst(x->x==ig_["TX27"],pbm.congrps)] = Float64(-1.411897512)
        end
        if haskey(ix_,"TY27")
            pb.x0[ix_["TY27"]] = Float64(-.1148065151)
        else
            pb.y0[findfirst(x->x==ig_["TY27"],pbm.congrps)] = Float64(-.1148065151)
        end
        if haskey(ix_,"TZ27")
            pb.x0[ix_["TZ27"]] = Float64(.44915582738)
        else
            pb.y0[findfirst(x->x==ig_["TZ27"],pbm.congrps)] = Float64(.44915582738)
        end
        if haskey(ix_,"KA27")
            pb.x0[ix_["KA27"]] = Float64(5.95875E-8)
        else
            pb.y0[findfirst(x->x==ig_["KA27"],pbm.congrps)] = Float64(5.95875E-8)
        end
        if haskey(ix_,"KB27")
            pb.x0[ix_["KB27"]] = Float64(-2.48391E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB27"],pbm.congrps)] = Float64(-2.48391E-13)
        end
        if haskey(ix_,"F27")
            pb.x0[ix_["F27"]] = Float64(407.03024568)
        else
            pb.y0[findfirst(x->x==ig_["F27"],pbm.congrps)] = Float64(407.03024568)
        end
        if haskey(ix_,"RX30")
            pb.x0[ix_["RX30"]] = Float64(.02082242153)
        else
            pb.y0[findfirst(x->x==ig_["RX30"],pbm.congrps)] = Float64(.02082242153)
        end
        if haskey(ix_,"RY30")
            pb.x0[ix_["RY30"]] = Float64(-1.238434791)
        else
            pb.y0[findfirst(x->x==ig_["RY30"],pbm.congrps)] = Float64(-1.238434791)
        end
        if haskey(ix_,"RZ30")
            pb.x0[ix_["RZ30"]] = Float64(.01389314763)
        else
            pb.y0[findfirst(x->x==ig_["RZ30"],pbm.congrps)] = Float64(.01389314763)
        end
        if haskey(ix_,"TX30")
            pb.x0[ix_["TX30"]] = Float64(-1.049686225)
        else
            pb.y0[findfirst(x->x==ig_["TX30"],pbm.congrps)] = Float64(-1.049686225)
        end
        if haskey(ix_,"TY30")
            pb.x0[ix_["TY30"]] = Float64(-.1299513286)
        else
            pb.y0[findfirst(x->x==ig_["TY30"],pbm.congrps)] = Float64(-.1299513286)
        end
        if haskey(ix_,"TZ30")
            pb.x0[ix_["TZ30"]] = Float64(.33798380231)
        else
            pb.y0[findfirst(x->x==ig_["TZ30"],pbm.congrps)] = Float64(.33798380231)
        end
        if haskey(ix_,"KA30")
            pb.x0[ix_["KA30"]] = Float64(4.5673127E-8)
        else
            pb.y0[findfirst(x->x==ig_["KA30"],pbm.congrps)] = Float64(4.5673127E-8)
        end
        if haskey(ix_,"KB30")
            pb.x0[ix_["KB30"]] = Float64(-1.79243E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB30"],pbm.congrps)] = Float64(-1.79243E-13)
        end
        if haskey(ix_,"F30")
            pb.x0[ix_["F30"]] = Float64(405.91764962)
        else
            pb.y0[findfirst(x->x==ig_["F30"],pbm.congrps)] = Float64(405.91764962)
        end
        if haskey(ix_,"RX37")
            pb.x0[ix_["RX37"]] = Float64(.01658816461)
        else
            pb.y0[findfirst(x->x==ig_["RX37"],pbm.congrps)] = Float64(.01658816461)
        end
        if haskey(ix_,"RY37")
            pb.x0[ix_["RY37"]] = Float64(-1.247226838)
        else
            pb.y0[findfirst(x->x==ig_["RY37"],pbm.congrps)] = Float64(-1.247226838)
        end
        if haskey(ix_,"RZ37")
            pb.x0[ix_["RZ37"]] = Float64(.01846788123)
        else
            pb.y0[findfirst(x->x==ig_["RZ37"],pbm.congrps)] = Float64(.01846788123)
        end
        if haskey(ix_,"TX37")
            pb.x0[ix_["TX37"]] = Float64(-.8617315756)
        else
            pb.y0[findfirst(x->x==ig_["TX37"],pbm.congrps)] = Float64(-.8617315756)
        end
        if haskey(ix_,"TY37")
            pb.x0[ix_["TY37"]] = Float64(-.1321089362)
        else
            pb.y0[findfirst(x->x==ig_["TY37"],pbm.congrps)] = Float64(-.1321089362)
        end
        if haskey(ix_,"TZ37")
            pb.x0[ix_["TZ37"]] = Float64(.28256800868)
        else
            pb.y0[findfirst(x->x==ig_["TZ37"],pbm.congrps)] = Float64(.28256800868)
        end
        if haskey(ix_,"KA37")
            pb.x0[ix_["KA37"]] = Float64(4.7465711E-8)
        else
            pb.y0[findfirst(x->x==ig_["KA37"],pbm.congrps)] = Float64(4.7465711E-8)
        end
        if haskey(ix_,"KB37")
            pb.x0[ix_["KB37"]] = Float64(-1.50881E-13)
        else
            pb.y0[findfirst(x->x==ig_["KB37"],pbm.congrps)] = Float64(-1.50881E-13)
        end
        if haskey(ix_,"F37")
            pb.x0[ix_["F37"]] = Float64(404.73590637)
        else
            pb.y0[findfirst(x->x==ig_["F37"],pbm.congrps)] = Float64(404.73590637)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eE", iet_)
        loaset(elftv,it,1,"RX")
        loaset(elftv,it,2,"RY")
        loaset(elftv,it,3,"RZ")
        loaset(elftv,it,4,"X")
        loaset(elftv,it,5,"Y")
        loaset(elftv,it,6,"Z")
        loaset(elftv,it,7,"TX")
        loaset(elftv,it,8,"TY")
        loaset(elftv,it,9,"TZ")
        loaset(elftv,it,10,"KA")
        loaset(elftv,it,11,"KB")
        loaset(elftv,it,12,"F")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"YRES")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EX1"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY1"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EX2"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY2"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F2"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EX3"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY3"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F4"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EX4"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY4"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F27"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EX5"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY5"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F30"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EX6"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        ename = "EY6"
        ie,ie_,_  = s2x_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eE")
        arrset(ielftype, ie, iet_["eE"])
        vname = "X1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Y1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "Z1"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RX37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RY37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "RZ37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="RZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TX37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TY37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "TZ37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="TZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KA37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KA",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "KB37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="KB",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "F37"
        iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="F",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="YRES",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["nuobservs"])
            ig = ig_["RX"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EX"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["RY"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EY"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "NOR2-MN-57-12"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
