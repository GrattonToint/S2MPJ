function HS114(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS114
#    *********
# 
#    An alkylation process problem.
# 
#    Source: problem 114 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: J.M. Collin, Jan 1990.
# 
#    classification = "C-QOR2-MY-10-11"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS114"

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
        v_["A"] = 0.99
        v_["B"] = 0.9
        v_["-A"] = -1.0*v_["A"]
        v_["-B"] = -1.0*v_["B"]
        v_["INVA"] = 1.0/v_["A"]
        v_["INVB"] = 1.0/v_["B"]
        v_["-INVA"] = -1.0*v_["INVA"]
        v_["-INVB"] = -1.0*v_["INVB"]
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
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.035)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(10.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(3.36)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(5.04)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-0.222)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(v_["-B"])
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(3.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(v_["-A"])
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C3")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(0.222)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(v_["INVB"])
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C4")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(v_["INVA"])
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.12)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(v_["-A"])
        ig,ig_,_ = s2mpj_ii("C6",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C6")
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.098)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.325)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(v_["-A"])
        ig,ig_,_ = s2mpj_ii("C7",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C7")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.12)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(v_["INVA"])
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-1.098)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-0.325)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(v_["INVA"])
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C9")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.22)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C10")
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C11")
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        pbm.gconst[ig_["C1"]] = Float64(-35.82)
        pbm.gconst[ig_["C2"]] = Float64(133.0)
        pbm.gconst[ig_["C3"]] = Float64(35.82)
        pbm.gconst[ig_["C4"]] = Float64(-133.0)
        pbm.gconst[ig_["C6"]] = Float64(-57.425)
        pbm.gconst[ig_["C8"]] = Float64(57.425)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 0.00001
        pb.xupper[ix_["X1"]] = 2000.0
        pb.xlower[ix_["X2"]] = 0.00001
        pb.xupper[ix_["X2"]] = 16000.0
        pb.xlower[ix_["X3"]] = 0.00001
        pb.xupper[ix_["X3"]] = 120.0
        pb.xlower[ix_["X4"]] = 0.00001
        pb.xupper[ix_["X4"]] = 5000.0
        pb.xlower[ix_["X5"]] = 0.00001
        pb.xupper[ix_["X5"]] = 2000.0
        pb.xlower[ix_["X6"]] = 85.0
        pb.xupper[ix_["X6"]] = 93.0
        pb.xlower[ix_["X7"]] = 90.0
        pb.xupper[ix_["X7"]] = 95.0
        pb.xlower[ix_["X8"]] = 3.0
        pb.xupper[ix_["X8"]] = 12.0
        pb.xlower[ix_["X9"]] = 1.2
        pb.xupper[ix_["X9"]] = 4.0
        pb.xlower[ix_["X10"]] = 145.0
        pb.xupper[ix_["X10"]] = 162.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(1745.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(1745.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(12000.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(12000.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(110.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(110.0)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(3048.0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(3048.0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(1974.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(1974.0)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(89.2)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(89.2)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(92.8)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(92.8)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(8.0)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(8.0)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(3.6)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(3.6)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(145.0)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(145.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eWSQ", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"W")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftp,it,1,"W")
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftp,it,1,"W")
        it,iet_,_ = s2mpj_ii( "eRAP1", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"S")
        loaset(elftv,it,3,"T")
        loaset(elftv,it,4,"U")
        loaset(elftp,it,1,"W")
        it,iet_,_ = s2mpj_ii( "eRAP2", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.063))
        ename = "CE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.13167))
        ename = "CE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.00667))
        ename = "CE3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eWSQ")
        arrset(ielftype,ie,iet_["eWSQ"])
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.038))
        ename = "CE4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.13167))
        ename = "CE5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD2")
        arrset(ielftype,ie,iet_["ePROD2"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00667))
        ename = "CE6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eWSQ")
        arrset(ielftype,ie,iet_["eWSQ"])
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.038))
        ename = "CE7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eRAP1")
        arrset(ielftype,ie,iet_["eRAP1"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="S",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="T",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="W",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(98000.0))
        ename = "CE8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eRAP2")
        arrset(ielftype,ie,iet_["eRAP2"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
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
        ig = ig_["C5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["C6"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C7"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["CE5"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["C8"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C10"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE7"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C11"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CE8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1768.80696
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
        pb.pbclass = "C-QOR2-MY-10-11"
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

    elseif action == "eWSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0*pbm.elpar[iel_][1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[2]
            g_[2] = pbm.elpar[iel_][1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]
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

    elseif action == "ePROD2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = pbm.elpar[iel_][1]*EV_[1]*EV_[2]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[2]^2
            g_[2] = 2.0*pbm.elpar[iel_][1]*EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 2.0*pbm.elpar[iel_][1]*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*pbm.elpar[iel_][1]*EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eRAP1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        WX = pbm.elpar[iel_][1]*EV_[1]
        DENOM = EV_[2]*EV_[3]+1000.0*EV_[4]
        f_   = WX/DENOM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]/DENOM
            g_[2] = -WX*EV_[3]/DENOM^2
            g_[3] = -WX*EV_[2]/DENOM^2
            g_[4] = -1000.0*WX/DENOM^2
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = -pbm.elpar[iel_][1]*EV_[3]/DENOM^2
                H_[2,1] = H_[1,2]
                H_[1,3] = -pbm.elpar[iel_][1]*EV_[2]/DENOM^2
                H_[3,1] = H_[1,3]
                H_[1,4] = -1000.0*pbm.elpar[iel_][1]/DENOM^2
                H_[4,1] = H_[1,4]
                H_[2,2] = 2.0*WX*EV_[3]*EV_[3]/DENOM^3
                H_[2,3] = 2.0*WX*EV_[2]*EV_[3]/DENOM^3-WX/DENOM^2
                H_[3,2] = H_[2,3]
                H_[2,4] = 2000.0*WX*EV_[3]/DENOM^3
                H_[4,2] = H_[2,4]
                H_[3,3] = 2.0*WX*EV_[2]*EV_[2]/DENOM^3
                H_[3,4] = 2000.0*WX*EV_[2]/DENOM^3
                H_[4,3] = H_[3,4]
                H_[4,4] = 2000000.0*WX/DENOM^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eRAP2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[2,3] = U_[2,3]+1
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[1]/IV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/IV_[2]
            g_[2] = -IV_[1]/IV_[2]^2
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0/IV_[2]^2
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*IV_[1]/IV_[2]^3
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

