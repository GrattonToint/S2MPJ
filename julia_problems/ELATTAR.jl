function ELATTAR(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ELATTAR
#    *********
# 
#    A nonlinear minmax problem in six variables.
# 
#    The problem is nonconvex and has several local minima.
# 
#    Source: 
#    R.A. El-Attar, M. Vidyasagar and S.R.K. Dutta,
#    "An algorithm for l_1-approximation",
#    SINUM 16, pp.70-86, 1979.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-CLOR2-AN-7-102"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ELATTAR"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling ELATTAR.")
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
        v_["6"] = 6
        v_["51"] = 51
        v_["T"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["51"])
            v_["T"*string(I)] = v_["T"]
            v_["T"] = 0.1+v_["T"]
            v_["ETI"] = exp(v_["T"*string(I)])
            v_["Y"*string(I)] = 0.5*v_["ETI"]
            v_["-2TI"] = -2.0*v_["T"*string(I)]
            v_["E-2TI"] = exp(v_["-2TI"])
            v_["Y"*string(I)] = v_["Y"*string(I)]-v_["E-2TI"]
            v_["-3TI"] = -3.0*v_["T"*string(I)]
            v_["E-3TI"] = exp(v_["-3TI"])
            v_["E-3TI/2"] = 0.5*v_["E-3TI"]
            v_["Y"*string(I)] = v_["Y"*string(I)]+v_["E-3TI/2"]
            v_["-3TI/2"] = 0.5*v_["-3TI"]
            v_["E-3TI/2"] = exp(v_["-3TI/2"])
            v_["7TI"] = 7.0*v_["T"*string(I)]
            v_["S7TI"] = sin(v_["7TI"])
            v_["TT"] = v_["E-3TI/2"]*v_["S7TI"]
            v_["TT"] = 1.5*v_["TT"]
            v_["Y"*string(I)] = v_["Y"*string(I)]+v_["TT"]
            v_["5TI"] = 5.0*v_["T"*string(I)]
            v_["-5TI/2"] = -0.5*v_["5TI"]
            v_["E-5TI/2"] = exp(v_["-5TI/2"])
            v_["S5TI"] = sin(v_["5TI"])
            v_["TT"] = v_["E-5TI/2"]*v_["S5TI"]
            v_["Y"*string(I)] = v_["Y"*string(I)]+v_["TT"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["6"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["U"])
        push!(valA,Float64(1.0))
        for I = Int64(v_["1"]):Int64(v_["51"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"F"*string(I))
            push!(irA,ig)
            push!(icA,ix_["U"])
            push!(valA,Float64(-1.0))
            ig,ig_,_ = s2mpj_ii("MF"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"MF"*string(I))
            push!(irA,ig)
            push!(icA,ix_["U"])
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
        for I = Int64(v_["1"]):Int64(v_["51"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
            v_["-Y"*string(I)] = -1.0*v_["Y"*string(I)]
            pbm.gconst[ig_["MF"*string(I)]] = Float64(v_["-Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(-2.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(-2.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(-2.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(-2.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(7.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(7.0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(-2.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(-2.0)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eET1", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        it,iet_,_ = s2mpj_ii( "eET2", iet_)
        loaset(elftv,it,1,"V5")
        loaset(elftv,it,2,"V6")
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["51"])
            ename = "EL1"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eET1")
            arrset(ielftype,ie,iet_["eET1"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
            ename = "EL2"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eET2")
            arrset(ielftype,ie,iet_["eET2"])
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["51"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EL1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EL2"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            ig = ig_["MF"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EL1"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EL2"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution       
# LO SOLTN               0.1427066255
# LO SOLTN               74.206179244
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLOR2-AN-7-102"
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

    elseif action == "eET1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = -EV_[2]*pbm.elpar[iel_][1]
        B = EV_[3]*pbm.elpar[iel_][1]+EV_[4]
        EA = exp(A)
        CB = cos(B)
        SB = sin(B)
        EACB = EA*CB
        EASB = EA*SB
        V1EACB = EV_[1]*EACB
        V1EASB = EV_[1]*EASB
        f_   = V1EACB
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EACB
            g_[2] = -pbm.elpar[iel_][1]*V1EACB
            g_[3] = -pbm.elpar[iel_][1]*V1EASB
            g_[4] = -V1EASB
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = -pbm.elpar[iel_][1]*EACB
                H_[2,1] = H_[1,2]
                H_[1,3] = -pbm.elpar[iel_][1]*EASB
                H_[3,1] = H_[1,3]
                H_[1,4] = -EASB
                H_[4,1] = H_[1,4]
                H_[2,2] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*V1EACB
                H_[2,3] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*V1EASB
                H_[3,2] = H_[2,3]
                H_[2,4] = pbm.elpar[iel_][1]*V1EASB
                H_[4,2] = H_[2,4]
                H_[3,3] = -pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*V1EACB
                H_[3,4] = -pbm.elpar[iel_][1]*V1EACB
                H_[4,3] = H_[3,4]
                H_[4,4] = -V1EACB
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eET2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = -EV_[2]*pbm.elpar[iel_][1]
        EA = exp(A)
        B = EV_[1]*EA
        f_   = B
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EA
            g_[2] = -pbm.elpar[iel_][1]*B
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -pbm.elpar[iel_][1]*EA
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*B
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

