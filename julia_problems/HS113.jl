function HS113(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS113
#    *********
# 
#    Source: problem 80 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: J-M Collin, Mar 1990.
# 
#    classification = "C-QQR2-AN-10-8"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS113"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 10
        v_["1"] = 1
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
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-14.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-16.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-20.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-40.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-6.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-154.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-40.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-14.0)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-4.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-5.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(3.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-9.0)
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-10.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(8.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(17.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-2.0)
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C3")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(8.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-5.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(2.0)
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C4")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(12.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(24.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(7.0)
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-8.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(12.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(2.0)
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C6")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(8.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(16.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C7")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(8.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-14.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(6.0)
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(3.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-6.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(192.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(7.0)
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
        pbm.gconst[ig_["OBJ"]] = Float64(-1352.0)
        pbm.gconst[ig_["C1"]] = Float64(-105.0)
        pbm.gconst[ig_["C3"]] = Float64(-12.0)
        pbm.gconst[ig_["C4"]] = Float64(-72.0)
        pbm.gconst[ig_["C5"]] = Float64(-4.0)
        pbm.gconst[ig_["C6"]] = Float64(34.0)
        pbm.gconst[ig_["C7"]] = Float64(8.0)
        pbm.gconst[ig_["C8"]] = Float64(768.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(5.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(5.0)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(5.0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(5.0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(7.0)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(7.0)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(10.0)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(10.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE9"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE10"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE11"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "CE12"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE2"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE4"])
        loaset(pbm.grelw,ig,posel,Float64(4.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE6"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE8"])
        loaset(pbm.grelw,ig,posel,Float64(7.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE9"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE10"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-3.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE2"])
        loaset(pbm.grelw,ig,posel,Float64(-4.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        ig = ig_["C5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-5.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE5"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C6"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE7"])
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-3.0))
        ig = ig_["C7"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE9"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE10"])
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE11"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        ig = ig_["C8"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE12"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-12.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               24.3062091
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QQR2-AN-10-8"
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

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

