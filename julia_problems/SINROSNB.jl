function SINROSNB(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SINROSNB
#    --------
# 
#    A variation on the extended Rosenbrock function in which
#    the squares are replaced by sines.
# 
#    Source: a modification of an original idea by
#    Ali Bouriacha, private communication.
# 
#    SIF input: Nick Gould and Ph. Toint, October, 1993.
# 
#    classification = "C-COQR2-AN-V-V"
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER     original value
# IE N                   4              $-PARAMETER
# IE N                   10             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SINROSNB"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling SINROSNB.")
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
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["PI"] = 3.1415926535
        v_["-PI"] = -1.0*v_["PI"]
        v_["-3PI/2"] = -1.5*v_["PI"]
        v_["2PI"] = 2.0*v_["PI"]
        v_["1-3PI/2"] = 1.0+v_["-3PI/2"]
        v_["1"] = 1
        v_["2"] = 2
        v_["N-1"] = -1+v_["N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("SQ1",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        for I = Int64(v_["2"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("SQ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
            arrset(pbm.gscale,ig,Float64(0.01))
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(I))
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
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
        pbm.gconst[ig_["SQ1"]] = Float64(v_["1-3PI/2"])
        for I = Int64(v_["2"]):Int64(v_["N"])
            pbm.gconst[ig_["SQ"*string(I)]] = Float64(v_["-3PI/2"])
            pbm.gconst[ig_["C"*string(I)]] = Float64(v_["-PI"])
        end
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for I = Int64(v_["2"]):Int64(v_["N"])
            arrset(grange,ig_["C"*string(I)],Float64(v_["2PI"]))
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = v_["-PI"]
        pb.xupper[ix_["X1"]] = v_["PI"]
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(10.0),pb.n)
        pb.x0[ix_["X1"]] = Float64(-1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eETYPE", iet_)
        loaset(elftv,it,1,"V1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "XSQ"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eETYPE")
                arrset(ielftype,ie,iet_["eETYPE"])
            end
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(10.0))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSIN",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["SQ1"]
        arrset(pbm.grftype,ig,"gSIN")
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig = ig_["SQ"*string(I)]
            arrset(pbm.grftype,ig,"gSIN")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["C"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution ( 1, 1, ..., 1 )
# LO SOLTN               0.0
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COQR2-AN-V-V"
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

    elseif action == "eETYPE"

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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSIN"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SING = sin(GVAR_)
        f_= SING+1.0
        if nargout>1
            g_ = cos(GVAR_)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -SING
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

