function ERRINROSNE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ERRINROSNE
#    --------
# 
#    A nonlinear function similar to the chained Rosenbrock
#    problem CHNROSNB.
#    This is a nonlinear equation variant of ERRINROS
# 
#    Source:
#    An error in specifying problem CHNROSNB.
#    SIF input: Ph. Toint, Sept 1990.
#               Nick Gould (nonlinear equation version), Jan 2019
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#    Number of variables (at most 50)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ERRINROSNE"

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
#       Alternative values for the SIF file parameters:
# IE N                   25             $-PARAMETER
# IE N                   50             $-PARAMETER
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
        v_["1"] = 1
        v_["2"] = 2
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
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("SQ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"SQ"*string(I))
            push!(irA,ig)
            push!(icA,ix_["X"*string(Int64(v_["I-1"]))])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("B"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"B"*string(I))
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(1.0))
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
        for I = Int64(v_["2"]):Int64(v_["N"])
            pbm.gconst[ig_["B"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(-1.0)
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)] = Float64(-1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eETYPE", iet_)
        loaset(elftv,it,1,"V1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["N"])
            ename = "ELA"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eETYPE")
                arrset(ielftype,ie,iet_["eETYPE"])
            end
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["AI2"] = v_["ALPH"*string(I)]*v_["ALPH"*string(I)]
            v_["AI"] = 16.0*v_["AI2"]
            ig = ig_["SQ"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ELA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["AI"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(10)            6.69463214
# LO SOLTN(25)            18.4609060
# LO SOLTN(50)            39.9041540
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-AN-V-V"
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
        f_   = -EV_[1]^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -2.0*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -2.0
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

