function LOTSCHD(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple quadratic program inspired by the economic lot scheduling
#    problem.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: T. Kuan, Virginia Tech., Spring 1993.
# 
#    classification = "C-CQLR2-AN-12-7"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LOTSCHD"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling LOTSCHD.")
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
        v_["X1"] = 1.502
        v_["X2"] = 1.126
        v_["X3"] = 0.815
        v_["X4"] = 1.268
        v_["X5"] = 1.502
        v_["X6"] = 0.740
        v_["A1"] = 1.8
        v_["A2"] = 3.2
        v_["A3"] = 6.1
        v_["A4"] = 3.2
        v_["A5"] = 1.8
        v_["A6"] = 7.4
        v_["C1"] = 11.0
        v_["C2"] = 3.0
        v_["C3"] = 20.0
        v_["C4"] = 17.0
        v_["C5"] = 9.0
        v_["C6"] = 20.0
        v_["C7"] = 126.1
        v_["1"] = 1
        v_["2"] = 2
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["6"])
            iv,ix_,_ = s2mpj_ii("T"*string(I),ix_)
            arrset(pb.xnames,iv,"T"*string(I))
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["6"])
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(v_["X"*string(I)]))
            ig,ig_,_ = s2mpj_ii("CONS7",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS7")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("CONS"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS"*string(I))
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(v_["A"*string(I)]))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(-1.0))
        end
        for I = Int64(v_["2"]):Int64(v_["4"])
            ig,ig_,_ = s2mpj_ii("CONS"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS"*string(I))
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(-1.0))
        end
        ig,ig_,_ = s2mpj_ii("CONS2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CONS2")
        push!(irA,ig)
        push!(icA,ix_["T3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["U3"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["2"])
            ig,ig_,_ = s2mpj_ii("CONS3",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS3")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(-1.0))
        end
        for I = Int64(v_["4"]):Int64(v_["6"])
            ig,ig_,_ = s2mpj_ii("CONS3",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS3")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(-1.0))
        end
        ig,ig_,_ = s2mpj_ii("CONS4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CONS4")
        push!(irA,ig)
        push!(icA,ix_["T1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["U1"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["5"]):Int64(v_["6"])
            ig,ig_,_ = s2mpj_ii("CONS4",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS4")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
            push!(valA,Float64(-1.0))
        end
        ig,ig_,_ = s2mpj_ii("CONS5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CONS5")
        push!(irA,ig)
        push!(icA,ix_["T1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["U1"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["5"])
            ig,ig_,_ = s2mpj_ii("CONS6",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CONS6")
            push!(irA,ig)
            push!(icA,ix_["T"*string(I)])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["U"*string(I)])
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
        for I = Int64(v_["1"]):Int64(v_["7"])
            pbm.gconst[ig_["CONS"*string(I)]] = Float64(v_["C"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["6"])
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gSQUARE")
        end
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-CQLR2-AN-12-7"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQUARE"

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

