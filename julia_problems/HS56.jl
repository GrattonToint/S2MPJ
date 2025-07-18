function HS56(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS56
#    *********
# 
#    Source: problem 56 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, April 1990
# 
#    classification = "C-COOR2-AN-7-4"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS56"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling HS56.")
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
        v_["N"] = 7
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["PAR"*string(Int64(v_["1"]))] = 4.2
        v_["PAR"*string(Int64(v_["2"]))] = 4.2
        v_["PAR"*string(Int64(v_["3"]))] = 4.2
        v_["PAR"*string(Int64(v_["4"]))] = 7.2
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
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["1"])))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["2"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["2"])))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["3"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["3"])))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["4"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["4"])))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(2.0))
        ig,ig_,_ = s2mpj_ii("CON"*string(Int64(v_["4"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON"*string(Int64(v_["4"])))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(2.0))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(0.50973968)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(0.50973968)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(0.50973968)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(0.50973968)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(0.50973968)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(0.50973968)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(0.98511078)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(0.98511078)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en3PROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        it,iet_,_ = s2mpj_ii( "ePSNSQ", iet_)
        loaset(elftv,it,1,"V1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en3PROD")
        arrset(ielftype,ie,iet_["en3PROD"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["1"]):Int64(v_["4"])
            v_["J"] = 1+I
            v_["K"] = 3+I
            ename = "E"*string(Int64(v_["J"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePSNSQ")
            arrset(ielftype,ie,iet_["ePSNSQ"])
            ename = "E"*string(Int64(v_["J"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["J"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="P",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["PAR"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["4"])
            v_["J"] = 1+I
            ig = ig_["CON"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["J"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                 -3.456
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
        pb.pbclass = "C-COOR2-AN-7-4"
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

    elseif action == "en3PROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]*EV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]*EV_[3]
            g_[2] = EV_[1]*EV_[3]
            g_[3] = EV_[1]*EV_[2]
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = EV_[3]
                H_[2,1] = H_[1,2]
                H_[1,3] = EV_[2]
                H_[3,1] = H_[1,3]
                H_[2,3] = EV_[1]
                H_[3,2] = H_[2,3]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePSNSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SUM = EV_[1]+EV_[1]
        f_   = -pbm.elpar[iel_][1]*sin(EV_[1])^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.elpar[iel_][1]*sin(SUM)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -2.0*pbm.elpar[iel_][1]*cos(SUM)
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

