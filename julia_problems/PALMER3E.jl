function PALMER3E(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER3E
#    *********
# 
#    A nonlinear least squares problem
#    arising from chemical kinetics.
# 
#     model: H-N=C=S TZVP + MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + L * EXP( -K X**2 )
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "C-CSBR2-RN-8-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER3E"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 23
        v_["1"] = 1
        v_["X1"] = -1.658063
        v_["X2"] = -1.570796
        v_["X3"] = -1.396263
        v_["X4"] = -1.221730
        v_["X5"] = -1.047198
        v_["X6"] = -0.872665
        v_["X7"] = -0.766531
        v_["X8"] = -0.698132
        v_["X9"] = -0.523599
        v_["X10"] = -0.349066
        v_["X11"] = -0.174533
        v_["X12"] = 0.0
        v_["X13"] = 0.174533
        v_["X14"] = 0.349066
        v_["X15"] = 0.523599
        v_["X16"] = 0.698132
        v_["X17"] = 0.766531
        v_["X18"] = 0.872665
        v_["X19"] = 1.047198
        v_["X20"] = 1.221730
        v_["X21"] = 1.396263
        v_["X22"] = 1.570796
        v_["X23"] = 1.658063
        v_["Y1"] = 64.87939
        v_["Y2"] = 50.46046
        v_["Y3"] = 28.2034
        v_["Y4"] = 13.4575
        v_["Y5"] = 4.6547
        v_["Y6"] = 0.59447
        v_["Y7"] = 0.0000
        v_["Y8"] = 0.2177
        v_["Y9"] = 2.3029
        v_["Y10"] = 5.5191
        v_["Y11"] = 8.5519
        v_["Y12"] = 9.8919
        v_["Y13"] = 8.5519
        v_["Y14"] = 5.5191
        v_["Y15"] = 2.3029
        v_["Y16"] = 0.2177
        v_["Y17"] = 0.0000
        v_["Y18"] = 0.59447
        v_["Y19"] = 4.6547
        v_["Y20"] = 13.4575
        v_["Y21"] = 28.2034
        v_["Y22"] = 50.46046
        v_["Y23"] = 64.87939
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("A0",ix_)
        arrset(pb.xnames,iv,"A0")
        iv,ix_,_ = s2mpj_ii("A2",ix_)
        arrset(pb.xnames,iv,"A2")
        iv,ix_,_ = s2mpj_ii("A4",ix_)
        arrset(pb.xnames,iv,"A4")
        iv,ix_,_ = s2mpj_ii("A6",ix_)
        arrset(pb.xnames,iv,"A6")
        iv,ix_,_ = s2mpj_ii("A8",ix_)
        arrset(pb.xnames,iv,"A8")
        iv,ix_,_ = s2mpj_ii("A10",ix_)
        arrset(pb.xnames,iv,"A10")
        iv,ix_,_ = s2mpj_ii("K",ix_)
        arrset(pb.xnames,iv,"K")
        iv,ix_,_ = s2mpj_ii("L",ix_)
        arrset(pb.xnames,iv,"L")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            v_["XQUART"] = v_["XSQR"]*v_["XSQR"]
            v_["X**6"] = v_["XSQR"]*v_["XQUART"]
            v_["X**8"] = v_["XSQR"]*v_["X**6"]
            v_["X**10"] = v_["XSQR"]*v_["X**8"]
            v_["X**12"] = v_["XSQR"]*v_["X**10"]
            v_["X**14"] = v_["XSQR"]*v_["X**12"]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
            push!(irA,ig)
            push!(icA,ix_["A0"])
            push!(valA,Float64(1.0))
            push!(irA,ig)
            push!(icA,ix_["A2"])
            push!(valA,Float64(v_["XSQR"]))
            push!(irA,ig)
            push!(icA,ix_["A4"])
            push!(valA,Float64(v_["XQUART"]))
            push!(irA,ig)
            push!(icA,ix_["A6"])
            push!(valA,Float64(v_["X**6"]))
            push!(irA,ig)
            push!(icA,ix_["A8"])
            push!(valA,Float64(v_["X**8"]))
            push!(irA,ig)
            push!(icA,ix_["A10"])
            push!(valA,Float64(v_["X**10"]))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["O"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["A0"]] = -Inf
        pb.xupper[ix_["A0"]] = +Inf
        pb.xlower[ix_["A2"]] = -Inf
        pb.xupper[ix_["A2"]] = +Inf
        pb.xlower[ix_["A4"]] = -Inf
        pb.xupper[ix_["A4"]] = +Inf
        pb.xlower[ix_["A6"]] = -Inf
        pb.xupper[ix_["A6"]] = +Inf
        pb.xlower[ix_["A8"]] = -Inf
        pb.xupper[ix_["A8"]] = +Inf
        pb.xlower[ix_["A10"]] = -Inf
        pb.xupper[ix_["A10"]] = +Inf
        pb.xlower[ix_["L"]] = -Inf
        pb.xupper[ix_["L"]] = +Inf
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"K")
        loaset(elftv,it,2,"L")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"XSQR")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "K"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="K",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="L",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="XSQR",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["XSQR"]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["O"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               5.0741053D-05
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSBR2-RN-8-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPON = exp(-EV_[1]*pbm.elpar[iel_][1])
        f_   = EV_[2]*EXPON
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.elpar[iel_][1]*EV_[2]*EXPON
            g_[2] = EXPON
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EV_[2]*EXPON
                H_[1,2] = -pbm.elpar[iel_][1]*EXPON
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

