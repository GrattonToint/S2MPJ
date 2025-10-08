function TOINTQOR(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TOINTQOR
#    *********
# 
#    Toint's  Quadratic Operations Research problem
# 
#    Source:
#    Ph.L. Toint,
#    "Some numerical results using a sparse matrix updating formula in
#    unconstrained optimization",
#    Mathematics of Computation 32(1):839-852, 1978.
# 
#    See also Buckley#55 (p.94) (With a slightly lower optimal value?)
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CQUR2-MN-50-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 8 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TOINTQOR"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling TOINTQOR.")
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
        v_["N"] = 50
        v_["ALPH1"] = 1.25
        v_["ALPH2"] = 1.40
        v_["ALPH3"] = 2.40
        v_["ALPH4"] = 1.40
        v_["ALPH5"] = 1.75
        v_["ALPH6"] = 1.20
        v_["ALPH7"] = 2.25
        v_["ALPH8"] = 1.20
        v_["ALPH9"] = 1.00
        v_["ALPH10"] = 1.10
        v_["ALPH11"] = 1.50
        v_["ALPH12"] = 1.60
        v_["ALPH13"] = 1.25
        v_["ALPH14"] = 1.25
        v_["ALPH15"] = 1.20
        v_["ALPH16"] = 1.20
        v_["ALPH17"] = 1.40
        v_["ALPH18"] = 0.50
        v_["ALPH19"] = 0.50
        v_["ALPH20"] = 1.25
        v_["ALPH21"] = 1.80
        v_["ALPH22"] = 0.75
        v_["ALPH23"] = 1.25
        v_["ALPH24"] = 1.40
        v_["ALPH25"] = 1.60
        v_["ALPH26"] = 2.00
        v_["ALPH27"] = 1.00
        v_["ALPH28"] = 1.60
        v_["ALPH29"] = 1.25
        v_["ALPH30"] = 2.75
        v_["ALPH31"] = 1.25
        v_["ALPH32"] = 1.25
        v_["ALPH33"] = 1.25
        v_["ALPH34"] = 3.00
        v_["ALPH35"] = 1.50
        v_["ALPH36"] = 2.00
        v_["ALPH37"] = 1.25
        v_["ALPH38"] = 1.40
        v_["ALPH39"] = 1.80
        v_["ALPH40"] = 1.50
        v_["ALPH41"] = 2.20
        v_["ALPH42"] = 1.40
        v_["ALPH43"] = 1.50
        v_["ALPH44"] = 1.25
        v_["ALPH45"] = 2.00
        v_["ALPH46"] = 1.50
        v_["ALPH47"] = 1.25
        v_["ALPH48"] = 1.40
        v_["ALPH49"] = 0.60
        v_["ALPH50"] = 1.50
        v_["BETA1"] = 1.0
        v_["BETA2"] = 1.5
        v_["BETA3"] = 1.0
        v_["BETA4"] = 0.1
        v_["BETA5"] = 1.5
        v_["BETA6"] = 2.0
        v_["BETA7"] = 1.0
        v_["BETA8"] = 1.5
        v_["BETA9"] = 3.0
        v_["BETA10"] = 2.0
        v_["BETA11"] = 1.0
        v_["BETA12"] = 3.0
        v_["BETA13"] = 0.1
        v_["BETA14"] = 1.5
        v_["BETA15"] = 0.15
        v_["BETA16"] = 2.0
        v_["BETA17"] = 1.0
        v_["BETA18"] = 0.1
        v_["BETA19"] = 3.0
        v_["BETA20"] = 0.1
        v_["BETA21"] = 1.2
        v_["BETA22"] = 1.0
        v_["BETA23"] = 0.1
        v_["BETA24"] = 2.0
        v_["BETA25"] = 1.2
        v_["BETA26"] = 3.0
        v_["BETA27"] = 1.5
        v_["BETA28"] = 3.0
        v_["BETA29"] = 2.0
        v_["BETA30"] = 1.0
        v_["BETA31"] = 1.2
        v_["BETA32"] = 2.0
        v_["BETA33"] = 1.0
        v_["D1"] = -5.0
        v_["D2"] = -5.0
        v_["D3"] = -5.0
        v_["D4"] = -2.5
        v_["D5"] = -6.0
        v_["D6"] = -6.0
        v_["D7"] = -5.0
        v_["D8"] = -6.0
        v_["D9"] = -10.0
        v_["D10"] = -6.0
        v_["D11"] = -5.0
        v_["D12"] = -9.0
        v_["D13"] = -2.0
        v_["D14"] = -7.0
        v_["D15"] = -2.5
        v_["D16"] = -6.0
        v_["D17"] = -5.0
        v_["D18"] = -2.0
        v_["D19"] = -9.0
        v_["D20"] = -2.0
        v_["D21"] = -5.0
        v_["D22"] = -5.0
        v_["D23"] = -2.5
        v_["D24"] = -5.0
        v_["D25"] = -6.0
        v_["D26"] = -10.0
        v_["D27"] = -7.0
        v_["D28"] = -10.0
        v_["D29"] = -6.0
        v_["D30"] = -5.0
        v_["D31"] = -4.0
        v_["D32"] = -4.0
        v_["D33"] = -4.0
        v_["1"] = 1
        v_["33"] = 33
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("GA"*string(I),ig_)
            arrset(gtype,ig,"<>")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(1.0))
            v_["SCALE"] = 1.0/v_["ALPH"*string(I)]
            arrset(pbm.gscale,ig,Float64(v_["SCALE"]))
        end
        ig,ig_,_ = s2mpj_ii("GB1",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X31"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB2",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB3",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB4",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB5",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB6",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB7",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB8",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X12"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB9",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X11"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X13"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X14"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X17"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB10",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X16"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X18"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB11",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X18"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB12",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X20"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X21"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("GB13",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X19"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB14",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X23"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X25"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X26"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB15",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X25"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X27"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X28"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB16",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X28"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X29"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X30"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB17",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X29"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X31"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X32"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB18",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X32"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X33"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X34"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB19",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X33"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X35"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB20",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X35"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X21"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X36"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB21",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X36"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X37"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X38"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB22",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X30"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X37"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X39"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB23",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X38"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X39"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X40"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB24",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X40"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X41"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X42"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB25",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X41"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X43"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X44"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X50"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB26",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X44"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X45"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X46"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X47"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB27",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X46"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X48"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB28",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X42"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X45"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X48"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X50"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X49"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("GB29",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X26"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X34"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X43"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("GB30",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X15"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X17"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X24"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["X47"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("GB31",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X49"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("GB32",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X22"])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("GB33",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X27"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["33"])
            v_["SCALE"] = 1.0/v_["BETA"*string(I)]
            ig,ig_,_ = s2mpj_ii("GB"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["SCALE"]))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["33"])
            pbm.gconst[ig_["GB"*string(I)]] = Float64(v_["D"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              1175.4722221
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CQUR2-MN-50-0"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

