function RES(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RES
#    *********
# 
#    Dassault France ressort (spring) problem
# 
#    SIF input:  A. R. Conn, June 1993.
# 
#    classification = "C-NLR2-MN-20-14"
# 
# 
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "RES"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("L0",ix_)
        arrset(pb.xnames,iv,"L0")
        iv,ix_,_ = s2mpj_ii("N",ix_)
        arrset(pb.xnames,iv,"N")
        iv,ix_,_ = s2mpj_ii("F",ix_)
        arrset(pb.xnames,iv,"F")
        iv,ix_,_ = s2mpj_ii("K",ix_)
        arrset(pb.xnames,iv,"K")
        iv,ix_,_ = s2mpj_ii("LB",ix_)
        arrset(pb.xnames,iv,"LB")
        iv,ix_,_ = s2mpj_ii("L",ix_)
        arrset(pb.xnames,iv,"L")
        iv,ix_,_ = s2mpj_ii("DE",ix_)
        arrset(pb.xnames,iv,"DE")
        iv,ix_,_ = s2mpj_ii("DI",ix_)
        arrset(pb.xnames,iv,"DI")
        iv,ix_,_ = s2mpj_ii("TO",ix_)
        arrset(pb.xnames,iv,"TO")
        iv,ix_,_ = s2mpj_ii("TOB",ix_)
        arrset(pb.xnames,iv,"TOB")
        iv,ix_,_ = s2mpj_ii("NU",ix_)
        arrset(pb.xnames,iv,"NU")
        iv,ix_,_ = s2mpj_ii("D",ix_)
        arrset(pb.xnames,iv,"D")
        iv,ix_,_ = s2mpj_ii("P",ix_)
        arrset(pb.xnames,iv,"P")
        iv,ix_,_ = s2mpj_ii("E",ix_)
        arrset(pb.xnames,iv,"E")
        iv,ix_,_ = s2mpj_ii("P0",ix_)
        arrset(pb.xnames,iv,"P0")
        iv,ix_,_ = s2mpj_ii("G",ix_)
        arrset(pb.xnames,iv,"G")
        iv,ix_,_ = s2mpj_ii("DM",ix_)
        arrset(pb.xnames,iv,"DM")
        iv,ix_,_ = s2mpj_ii("FR",ix_)
        arrset(pb.xnames,iv,"FR")
        iv,ix_,_ = s2mpj_ii("TOLIM",ix_)
        arrset(pb.xnames,iv,"TOLIM")
        iv,ix_,_ = s2mpj_ii("TOBLIM",ix_)
        arrset(pb.xnames,iv,"TOBLIM")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("E1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E1")
        iv = ix_["F"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E2")
        iv = ix_["K"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E3")
        iv = ix_["DE"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["D"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["DM"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E4")
        iv = ix_["DI"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["D"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["DM"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E5")
        iv = ix_["D"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["P"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["E"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E6")
        iv = ix_["NU"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["N"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E7")
        iv = ix_["D"]
        pbm.A[ig,iv] += Float64(1.5)
        iv = ix_["L0"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E8",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E8")
        iv = ix_["L"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["LB"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["FR"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E9",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E9")
        iv = ix_["LB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E10")
        iv = ix_["L"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["L0"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["F"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E11",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E11")
        iv = ix_["TO"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E12",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"E12")
        iv = ix_["TOB"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("E13",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E13")
        iv = ix_["TO"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TOLIM"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("E14",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E14")
        iv = ix_["TOB"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["TOBLIM"]
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
        pbm.gconst[ig_["E6"]] = Float64(-2.0)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["L0"]] = 100.0
        pb.xupper[ix_["N"]] = 100.0
        pb.xupper[ix_["F"]] = 30.0
        pb.xupper[ix_["K"]] = 100.0
        pb.xupper[ix_["LB"]] = 50.0
        pb.xupper[ix_["L"]] = 50.0
        pb.xupper[ix_["DE"]] = 30.0
        pb.xupper[ix_["DI"]] = 30.0
        pb.xupper[ix_["TO"]] = 800.0
        pb.xupper[ix_["TOB"]] = 800.0
        pb.xupper[ix_["NU"]] = 50.0
        pb.xlower[ix_["NU"]] = 0.5
        pb.xupper[ix_["D"]] = 10.0
        pb.xlower[ix_["D"]] = 0.1
        pb.xupper[ix_["P"]] = 20.0
        pb.xupper[ix_["E"]] = 10.0
        pb.xupper[ix_["P0"]] = 1000.0
        pb.xlower[ix_["P0"]] = 1.0
        pb.xupper[ix_["G"]] = 80000.0
        pb.xlower[ix_["G"]] = 40000.0
        pb.xupper[ix_["DM"]] = 30.0
        pb.xlower[ix_["DM"]] = 0.1
        pb.xupper[ix_["FR"]] = 50.0
        pb.xupper[ix_["TOLIM"]] = 1000.0
        pb.xlower[ix_["TOLIM"]] = 100.0
        pb.xupper[ix_["TOBLIM"]] = 1000.0
        pb.xlower[ix_["TOBLIM"]] = 100.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["L0"]] = Float64(1.5000e-01)
        pb.x0[ix_["N"]] = Float64(2.4079e+01)
        pb.x0[ix_["F"]] = Float64(9.2459e-15)
        pb.x0[ix_["K"]] = Float64(0.0000)
        pb.x0[ix_["LB"]] = Float64(0.0000)
        pb.x0[ix_["L"]] = Float64(1.5000e-01)
        pb.x0[ix_["DE"]] = Float64(6.8120)
        pb.x0[ix_["DI"]] = Float64(6.6120)
        pb.x0[ix_["TO"]] = Float64(0.0000)
        pb.x0[ix_["TOB"]] = Float64(0.0000)
        pb.x0[ix_["NU"]] = Float64(2.2079e+01)
        pb.x0[ix_["D"]] = Float64(1.0000e-01)
        pb.x0[ix_["P"]] = Float64(6.5268e-01)
        pb.x0[ix_["E"]] = Float64(5.5268e-01)
        pb.x0[ix_["P0"]] = Float64(6.5887e+02)
        pb.x0[ix_["G"]] = Float64(6.5887e+04)
        pb.x0[ix_["DM"]] = Float64(6.7120)
        pb.x0[ix_["FR"]] = Float64(1.5000e-01)
        pb.x0[ix_["TOLIM"]] = Float64(1.0000e+02)
        pb.x0[ix_["TOBLIM"]] = Float64(1.0000e+02)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "en311d14", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        loaset(elftv,it,3,"X")
        loaset(elftv,it,4,"Y")
        loaset(elftv,it,5,"Z")
        it,iet_,_ = s2mpj_ii( "en14d31", iet_)
        loaset(elftv,it,1,"W")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"Y")
        loaset(elftv,it,4,"Z")
        it,iet_,_ = s2mpj_ii( "en11d3", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "en111d2", iet_)
        loaset(elftv,it,1,"W")
        loaset(elftv,it,2,"X")
        loaset(elftv,it,3,"Y")
        loaset(elftv,it,4,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EL1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en311d14")
        arrset(ielftype,ie,iet_["en311d14"])
        vname = "DM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NU"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "P0"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "G"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en14d31")
        arrset(ielftype,ie,iet_["en14d31"])
        vname = "G"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "DM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "NU"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "NU"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "P"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en11d3")
        arrset(ielftype,ie,iet_["en11d3"])
        vname = "P0"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "DM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EL5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en111d2")
        arrset(ielftype,ie,iet_["en111d2"])
        vname = "G"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "D"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "E"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "DM"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-NLR2-MN-20-14"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,3.1415926535)
        return pbm

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

    elseif action == "en311d14"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V3WX = EV_[1]^3*EV_[2]*EV_[3]
        YZ4 = EV_[4]*EV_[5]^4
        f_   = V3WX/YZ4
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (3.0*EV_[1]^2*EV_[2]*EV_[3])/YZ4
            g_[2] = (EV_[1]^3*EV_[3])/YZ4
            g_[3] = (EV_[1]^3*EV_[2])/YZ4
            g_[4] = -V3WX/(EV_[4]*YZ4)
            g_[5] = -(4.0*V3WX)/(EV_[5]*YZ4)
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,1] = (6.0*EV_[1]*EV_[2]*EV_[3])/YZ4
                H_[1,2] = (3.0*EV_[1]^2*EV_[3])/YZ4
                H_[2,1] = H_[1,2]
                H_[1,3] = (3.0*EV_[1]^2*EV_[2])/YZ4
                H_[3,1] = H_[1,3]
                H_[1,4] = -(3.0*EV_[1]^2*EV_[2])/(YZ4*EV_[4])
                H_[4,1] = H_[1,4]
                H_[1,5] = -(12.0*EV_[1]^2*EV_[2]*EV_[3])/(YZ4*EV_[5])
                H_[5,1] = H_[1,5]
                H_[2,3] = EV_[1]^3/YZ4
                H_[3,2] = H_[2,3]
                H_[2,4] = -(EV_[1]^3*EV_[3])/(YZ4*EV_[4])
                H_[4,2] = H_[2,4]
                H_[2,5] = -(4.0*EV_[1]^3*EV_[3])/(YZ4*EV_[5])
                H_[5,2] = H_[2,5]
                H_[3,4] = -(EV_[1]^3*EV_[2])/(YZ4*EV_[4])
                H_[4,3] = H_[3,4]
                H_[3,5] = -(4.0*EV_[1]^3*EV_[2])/(YZ4*EV_[5])
                H_[5,3] = H_[3,5]
                H_[4,4] = -(2.0*V3WX)/(EV_[4]^2*YZ4)
                H_[4,5] = (4.0*V3WX)/(EV_[4]*EV_[5]*YZ4)
                H_[5,4] = H_[4,5]
                H_[5,5] = (20.0*V3WX)/(EV_[5]^2*YZ4)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en14d31"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        WX4 = EV_[1]*EV_[2]^4
        Y3Z = EV_[3]^3*EV_[4]
        f_   = WX4/Y3Z
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]^4/Y3Z
            g_[2] = (4.0*EV_[1]*EV_[2]^3)/Y3Z
            g_[3] = -(3.0*WX4)/(EV_[3]*Y3Z)
            g_[4] = -WX4/(EV_[4]*Y3Z)
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = (4.0*EV_[2]^3)/Y3Z
                H_[2,1] = H_[1,2]
                H_[1,3] = -(3.0*EV_[2]^4)/(EV_[3]*Y3Z)
                H_[3,1] = H_[1,3]
                H_[1,4] = -EV_[2]^4/(EV_[4]*Y3Z)
                H_[4,1] = H_[1,4]
                H_[2,2] = (12.0*EV_[1]*EV_[2]^2)/Y3Z
                H_[2,3] = -(12.0*EV_[1]*EV_[2]^3)/(EV_[3]*Y3Z)
                H_[3,2] = H_[2,3]
                H_[2,4] = -(4.0*EV_[1]*EV_[2]^3)/(EV_[4]*Y3Z)
                H_[4,2] = H_[2,4]
                H_[3,3] = (12.0*WX4)/(EV_[3]^2*Y3Z)
                H_[3,4] = (3.0*WX4)/(EV_[4]*EV_[3]*Y3Z)
                H_[4,3] = H_[3,4]
                H_[4,4] = (2.0*WX4)/(EV_[4]^2*Y3Z)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en11d3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]*EV_[2])/(pbm.efpar[1]*EV_[3]^3)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]/(pbm.efpar[1]*EV_[3]^3)
            g_[2] = EV_[1]/(pbm.efpar[1]*EV_[3]^3)
            g_[3] = -(3.0*EV_[1]*EV_[2])/(pbm.efpar[1]*EV_[3]^4)
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = 1.0/(pbm.efpar[1]*EV_[3]^3)
                H_[2,1] = H_[1,2]
                H_[1,3] = -(3.0*EV_[2])/(pbm.efpar[1]*EV_[3]^4)
                H_[3,1] = H_[1,3]
                H_[2,3] = -(3.0*EV_[1])/(pbm.efpar[1]*EV_[3]^4)
                H_[3,2] = H_[2,3]
                H_[3,3] = (12.0*EV_[1]*EV_[2])/(pbm.efpar[1]*EV_[3]^5)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en111d2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]*EV_[2]*EV_[3])/(pbm.efpar[1]*EV_[4]^2)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (EV_[2]*EV_[3])/(pbm.efpar[1]*EV_[4]^2)
            g_[2] = (EV_[1]*EV_[3])/(pbm.efpar[1]*EV_[4]^2)
            g_[3] = (EV_[1]*EV_[2])/(pbm.efpar[1]*EV_[4]^2)
            g_[4] = -(2.0*EV_[1]*EV_[2]*EV_[3])/(pbm.efpar[1]*EV_[4]^3)
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = EV_[3]/(pbm.efpar[1]*EV_[4]^2)
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]/(pbm.efpar[1]*EV_[4]^2)
                H_[3,1] = H_[1,3]
                H_[1,4] = -(2.0*EV_[2]*EV_[3])/(pbm.efpar[1]*EV_[4]^3)
                H_[4,1] = H_[1,4]
                H_[2,3] = EV_[1]/(pbm.efpar[1]*EV_[4]^2)
                H_[3,2] = H_[2,3]
                H_[2,4] = -(2.0*EV_[1]*EV_[3])/(pbm.efpar[1]*EV_[4]^3)
                H_[4,2] = H_[2,4]
                H_[3,4] = -(2.0*EV_[1]*EV_[2])/(pbm.efpar[1]*EV_[4]^3)
                H_[4,3] = H_[3,4]
                H_[4,4] = (6.0*EV_[1]*EV_[2]*EV_[3])/(pbm.efpar[1]*EV_[4]^4)
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

