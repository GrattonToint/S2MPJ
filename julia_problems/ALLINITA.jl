function ALLINITA(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ALLINITA
#    *********
# 
#    A problem with "all in it". Intended to verify that changes
#    to LANCELOT are safe. Multiple constrained version
# 
#    Source:
#    N. Gould: private communication.
# 
#    SIF input: Nick Gould, March 2013.
# 
#    classification = "C-COOR2-AY-4-4"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ALLINITA"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling ALLINITA.")
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
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("FT1",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FT2",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("FT3",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FT4",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FT5",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("FT6",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FNT1",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FNT2",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("FNT3",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FNT4",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("FNT5",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("FNT6",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C2")
        ig,ig_,_ = s2mpj_ii("L1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L1")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("L2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"L2")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
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
        pbm.gconst[ig_["FT2"]] = Float64(1.0)
        pbm.gconst[ig_["FT5"]] = Float64(3.0)
        pbm.gconst[ig_["FNT2"]] = Float64(1.0)
        pbm.gconst[ig_["FNT5"]] = Float64(4.0)
        pbm.gconst[ig_["C1"]] = Float64(1.0)
        pbm.gconst[ig_["L1"]] = Float64(1.5)
        pbm.gconst[ig_["L2"]] = Float64(0.25)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = -Inf
        pb.xupper[ix_["X1"]] = +Inf
        pb.xlower[ix_["X2"]] = 1.0
        pb.xlower[ix_["X3"]] = -1.0e+10
        pb.xupper[ix_["X3"]] = 1.0
        pb.xlower[ix_["X4"]] = 2.0
        pb.xupper[ix_["X4"]] = 2.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eSQR2", iet_)
        loaset(elftv,it,1,"Y")
        loaset(elftv,it,2,"Z")
        it,iet_,_ = s2mpj_ii( "eSINSQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "ePRODSQR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "FT3E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FT4E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FT5E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FT4E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR2")
        arrset(ielftype,ie,iet_["eSQR2"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FT56E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINSQR")
        arrset(ielftype,ie,iet_["eSINSQR"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FT5E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODSQR")
        arrset(ielftype,ie,iet_["ePRODSQR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FNT3E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FNT4E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
        end
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FNT4E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR2")
        arrset(ielftype,ie,iet_["eSQR2"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FNT56E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSINSQR")
        arrset(ielftype,ie,iet_["eSINSQR"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "FNT5E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePRODSQR")
        arrset(ielftype,ie,iet_["ePRODSQR"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["FT3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT3E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["FT4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT4E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FT4E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["FT5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT56E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FT5E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["FT6"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT56E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["FNT1"]
        arrset(pbm.grftype,ig,"gL2")
        ig = ig_["FNT2"]
        arrset(pbm.grftype,ig,"gL2")
        ig = ig_["FNT3"]
        arrset(pbm.grftype,ig,"gL2")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FNT3E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["FNT4"]
        arrset(pbm.grftype,ig,"gL2")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FNT4E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FNT4E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["FNT5"]
        arrset(pbm.grftype,ig,"gL2")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FNT56E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FNT5E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["FNT6"]
        arrset(pbm.grftype,ig,"gL2")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FNT56E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT3E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FT4E1"])
        loaset(pbm.grelw,ig,posel, 1.)
        ig = ig_["C2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["FT4E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["FT5E1"])
        loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COOR2-AY-4-4"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQR"

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

    elseif action == "eSQR2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
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

    elseif action == "eSINSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SINX = sin(EV_[1])
        COSX = cos(EV_[1])
        f_   = SINX*SINX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*SINX*COSX
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0*(COSX*COSX-SINX*SINX)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePRODSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        XX = EV_[1]*EV_[1]
        YY = EV_[2]*EV_[2]
        f_   = XX*YY
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]*YY
            g_[2] = 2.0*XX*EV_[2]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*YY
                H_[1,2] = 4.0*EV_[1]*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*XX
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

