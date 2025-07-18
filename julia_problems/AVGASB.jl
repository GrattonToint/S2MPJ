function AVGASB(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : AVGASB
#    *********
# 
#    LP problem avgas, with pd Hessian and integer variables
# 
#    SIF input: S. Leyffer, May 1998.
# 
#    classification = "C-CQLR2-AN-8-10"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "AVGASB"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling AVGASB.")
    end

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
        v_["8"] = 8
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["8"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-2.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-3.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-2.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-4.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-3.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-5.0))
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON1")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON2")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON3")
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON4",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON4")
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON5")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON6",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON6")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON7",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON7")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON8",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON8")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(5.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(3.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-3.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("CON9",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON9")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-3.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-5.0))
        ig,ig_,_ = s2mpj_ii("CON10",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON10")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-3.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-2.0))
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
        pbm.gconst[ig_["CON1"]] = Float64(-1.0)
        pbm.gconst[ig_["CON2"]] = Float64(-1.0)
        pbm.gconst[ig_["CON3"]] = Float64(-1.0)
        pbm.gconst[ig_["CON4"]] = Float64(-1.0)
        pbm.gconst[ig_["CON5"]] = Float64(-2.0)
        pbm.gconst[ig_["CON6"]] = Float64(-2.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = fill(1.0,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V1")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "X1S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X2S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X3S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X4S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X5S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X6S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X7S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X8S"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X1X2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X2X3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X3X4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X4X5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X5X6"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X6X7"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X7X8"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X1S"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X2S"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3S"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4S"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X5S"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X6S"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X7S"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X8S"])
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X1X2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X2X3"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3X4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4X5"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X5X6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X6X7"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X7X8"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CQLR2-AN-8-10"
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]
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

