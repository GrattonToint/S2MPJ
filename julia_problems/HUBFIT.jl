function HUBFIT(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HUBFIT
#    *********
#    Variable dimension full rank linear problem
#    An elementary fit using the Huber loss function
# 
#    Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "The LANCELOT User's Manual",
#    Dept of Maths, FUNDP, 1991.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "C-COLR2-AN-2-1"
# 
#    Data points
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HUBFIT"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling HUBFIT.")
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
        v_["X1"] = 0.1
        v_["X2"] = 0.3
        v_["X3"] = 0.5
        v_["X4"] = 0.7
        v_["X5"] = 0.9
        v_["Y1"] = 0.25
        v_["Y2"] = 0.3
        v_["Y3"] = 0.625
        v_["Y4"] = 0.701
        v_["Y5"] = 1.0
        v_["C"] = 0.85
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("a",ix_)
        arrset(pb.xnames,iv,"a")
        iv,ix_,_ = s2mpj_ii("b",ix_)
        arrset(pb.xnames,iv,"b")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("Obj1",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(v_["X1"]))
        push!(irA,ig)
        push!(icA,ix_["b"])
        push!(valA,Float64(1.0))
        arrset(pbm.gscale,ig,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Obj2",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(v_["X2"]))
        push!(irA,ig)
        push!(icA,ix_["b"])
        push!(valA,Float64(1.0))
        arrset(pbm.gscale,ig,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Obj3",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(v_["X3"]))
        push!(irA,ig)
        push!(icA,ix_["b"])
        push!(valA,Float64(1.0))
        arrset(pbm.gscale,ig,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Obj4",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(v_["X4"]))
        push!(irA,ig)
        push!(icA,ix_["b"])
        push!(valA,Float64(1.0))
        arrset(pbm.gscale,ig,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Obj5",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(v_["X5"]))
        push!(irA,ig)
        push!(icA,ix_["b"])
        push!(valA,Float64(1.0))
        arrset(pbm.gscale,ig,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("Cons",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"Cons")
        push!(irA,ig)
        push!(icA,ix_["a"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["b"])
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
        pbm.gconst[ig_["Obj1"]] = Float64(v_["Y1"])
        pbm.gconst[ig_["Obj2"]] = Float64(v_["Y2"])
        pbm.gconst[ig_["Obj3"]] = Float64(v_["Y3"])
        pbm.gconst[ig_["Obj4"]] = Float64(v_["Y4"])
        pbm.gconst[ig_["Obj5"]] = Float64(v_["Y5"])
        pbm.gconst[ig_["Cons"]] = Float64(v_["C"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["b"]] = -Inf
        pb.xupper[ix_["b"]] = +Inf
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gHUBER",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["Obj1"]
        arrset(pbm.grftype,ig,"gHUBER")
        ig = ig_["Obj2"]
        arrset(pbm.grftype,ig,"gHUBER")
        ig = ig_["Obj3"]
        arrset(pbm.grftype,ig,"gHUBER")
        ig = ig_["Obj4"]
        arrset(pbm.grftype,ig,"gHUBER")
        ig = ig_["Obj5"]
        arrset(pbm.grftype,ig,"gHUBER")
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-COLR2-AN-2-1"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gHUBER"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        HUBERK = 1.5
        ABSA = abs(GVAR_)
        OUT = ABSA>HUBERK
        if OUT
            FF = HUBERK*ABSA-0.5*HUBERK*HUBERK
        end
        NEGOUT = OUT&&(GVAR_<0.0)
        POSOUT = OUT&&(GVAR_>=0.0)
        if POSOUT
            GG = HUBERK
        end
        if NEGOUT
            GG = -HUBERK
        end
        if OUT
            HH = 0.0
        end
        if !OUT
            FF = 0.5*ABSA*ABSA
        end
        if !OUT
            GG = GVAR_
        end
        if !OUT
            HH = 1.0
        end
        f_= FF
        if nargout>1
            g_ = GG
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = HH
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

