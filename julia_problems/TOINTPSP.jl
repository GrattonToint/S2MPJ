function TOINTPSP(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TOINTPSP
#    *********
# 
#    Toint's PSP Operations Research problem
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
#    classification = "C-OUR2-AN-50-0"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TOINTPSP"

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
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("GA"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            v_["SCALE"] = 1.0/v_["ALPH"*string(I)]
            arrset(pbm.gscale,ig,Float64(v_["SCALE"]))
        end
        ig,ig_,_ = s2mpj_ii("GB1",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X31"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB2",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB3",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB4",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB5",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB6",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB7",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB8",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB9",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB10",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB11",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB12",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GB13",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB14",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X25"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X26"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB15",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X25"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X27"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X28"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB16",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X28"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X29"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X30"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB17",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X29"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X31"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X32"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB18",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X32"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X33"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X34"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB19",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X33"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X35"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB20",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X35"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X36"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB21",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X36"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X37"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X38"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB22",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X30"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X37"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X39"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB23",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X38"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X39"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X40"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB24",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X40"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X41"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB25",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X41"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X50"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB26",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X44"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X47"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB27",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X46"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X48"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB28",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X42"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X45"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X48"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X50"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X49"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("GB29",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X26"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X34"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X43"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GB30",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X47"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GB31",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X49"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GB32",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("GB33",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X27"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gACT",igt_)
        it,igt_,_ = s2mpj_ii("gBBT",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["GA"*string(I)]
            arrset(pbm.grftype,ig,"gACT")
        end
        for I = Int64(v_["1"]):Int64(v_["33"])
            ig = ig_["GB"*string(I)]
            arrset(pbm.grftype,ig,"gBBT")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN              225.56040942
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-50-0"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gACT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= (GVAR_-5.0)^2
        if nargout>1
            g_ = 2.0*GVAR_-10.0
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

    elseif action == "gBBT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TPOS = GVAR_>=0.1
        if TPOS
            FF = 1.0/GVAR_
        end
        if !TPOS
            FF = 20.0-100.0*GVAR_
        end
        if TPOS
            GG = -1.0/GVAR_^2
        end
        if !TPOS
            GG = -100.0
        end
        if TPOS
            HH = 2.0/GVAR_^3
        end
        if !TPOS
            HH = 0.0
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

