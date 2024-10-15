function DECONVU(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DECONVU
#    *********
# 
#    A problem arising in deconvolution analysis 
#    (unconstrained version).
# 
#    Source:  
#    J.P. Rasson, Private communication, 1996.
# 
#    SIF input: Ph. Toint, Nov 1996.
#    unititialized variables fixed at zero, Nick Gould, Feb, 2013
# 
#    classification = "C-SXR2-MN-61-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DECONVU"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["0"] = 0
        v_["1"] = 1
        v_["LGSG"] = 11
        v_["LGTR"] = 40
        v_["-LGSG"] = -1*v_["LGSG"]
        v_["TR1"] = 0.0000000000
        v_["TR2"] = 0.0000000000
        v_["TR3"] = 1.600000E-03
        v_["TR4"] = 5.400000E-03
        v_["TR5"] = 7.020000E-02
        v_["TR6"] = 0.1876000000
        v_["TR7"] = 0.3320000000
        v_["TR8"] = 0.7640000000
        v_["TR9"] = 0.9320000000
        v_["TR10"] = 0.8120000000
        v_["TR11"] = 0.3464000000
        v_["TR12"] = 0.2064000000
        v_["TR13"] = 8.300000E-02
        v_["TR14"] = 3.400000E-02
        v_["TR15"] = 6.179999E-02
        v_["TR16"] = 1.2000000000
        v_["TR17"] = 1.8000000000
        v_["TR18"] = 2.4000000000
        v_["TR19"] = 9.0000000000
        v_["TR20"] = 2.4000000000
        v_["TR21"] = 1.8010000000
        v_["TR22"] = 1.3250000000
        v_["TR23"] = 7.620000E-02
        v_["TR24"] = 0.2104000000
        v_["TR25"] = 0.2680000000
        v_["TR26"] = 0.5520000000
        v_["TR27"] = 0.9960000000
        v_["TR28"] = 0.3600000000
        v_["TR29"] = 0.2400000000
        v_["TR30"] = 0.1510000000
        v_["TR31"] = 2.480000E-02
        v_["TR32"] = 0.2432000000
        v_["TR33"] = 0.3602000000
        v_["TR34"] = 0.4800000000
        v_["TR35"] = 1.8000000000
        v_["TR36"] = 0.4800000000
        v_["TR37"] = 0.3600000000
        v_["TR38"] = 0.2640000000
        v_["TR39"] = 6.000000E-03
        v_["TR40"] = 6.000000E-03
        v_["SSG1"] = 1.000000E-02
        v_["SSG2"] = 2.000000E-02
        v_["SSG3"] = 0.4000000000
        v_["SSG4"] = 0.6000000000
        v_["SSG5"] = 0.8000000000
        v_["SSG6"] = 3.0000000000
        v_["SSG7"] = 0.8000000000
        v_["SSG8"] = 0.6000000000
        v_["SSG9"] = 0.4400000000
        v_["SSG10"] = 1.000000E-02
        v_["SSG11"] = 1.000000E-02
        v_["CC1"] = 0.0
        v_["CC2"] = 0.0
        v_["CC3"] = 0.0
        v_["CC4"] = 0.0
        v_["CC5"] = 0.0
        v_["CC6"] = 0.0
        v_["CC7"] = 0.0
        v_["CC8"] = 0.0
        v_["CC9"] = 0.0
        v_["CC10"] = 0.0
        v_["CC11"] = 0.0
        v_["CC12"] = 0.0
        v_["CC13"] = 0.0
        v_["CC14"] = 0.0
        v_["CC15"] = 0.0
        v_["CC16"] = 0.0
        v_["CC17"] = 0.0
        v_["CC18"] = 0.0
        v_["CC19"] = 0.0
        v_["CC20"] = 0.0
        v_["CC21"] = 0.0
        v_["CC22"] = 0.0
        v_["CC23"] = 0.0
        v_["CC24"] = 0.0
        v_["CC25"] = 0.0
        v_["CC26"] = 0.0
        v_["CC27"] = 0.0
        v_["CC28"] = 0.0
        v_["CC29"] = 0.0
        v_["CC30"] = 0.0
        v_["CC31"] = 0.0
        v_["CC32"] = 0.0
        v_["CC33"] = 0.0
        v_["CC34"] = 0.0
        v_["CC35"] = 0.0
        v_["CC36"] = 0.0
        v_["CC37"] = 0.0
        v_["CC38"] = 0.0
        v_["CC39"] = 0.0
        v_["CC40"] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for K = Int64(v_["-LGSG"]):Int64(v_["LGTR"])
            iv,ix_,_ = s2mpj_ii("C"*string(K),ix_)
            arrset(pb.xnames,iv,"C"*string(K))
        end
        for I = Int64(v_["1"]):Int64(v_["LGSG"])
            iv,ix_,_ = s2mpj_ii("SG"*string(I),ix_)
            arrset(pb.xnames,iv,"SG"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for K = Int64(v_["1"]):Int64(v_["LGTR"])
            ig,ig_,_ = s2mpj_ii("R"*string(K),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for K = Int64(v_["1"]):Int64(v_["LGTR"])
            pbm.gconst[ig_["R"*string(K)]] = Float64(v_["TR"*string(K)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for K = Int64(v_["-LGSG"]):Int64(v_["0"])
            pb.xlower[ix_["C"*string(K)]] = 0.0
            pb.xupper[ix_["C"*string(K)]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for K = Int64(v_["1"]):Int64(v_["LGTR"])
            pb.x0[ix_["C"*string(K)]] = Float64(v_["CC"*string(K)])
        end
        for I = Int64(v_["1"]):Int64(v_["LGSG"])
            pb.x0[ix_["SG"*string(I)]] = Float64(v_["SSG"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"IDX")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["1"]):Int64(v_["LGTR"])
            for I = Int64(v_["1"]):Int64(v_["LGSG"])
                v_["K-I"] = K-I
                v_["K-I+1"] = 1+v_["K-I"]
                v_["RIDX"] = Float64(v_["K-I+1"])
                ename = "PROD"*string(K)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePR")
                arrset(ielftype,ie,iet_["ePR"])
                vname = "SG"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "C"*string(Int64(v_["K-I+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="IDX",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RIDX"]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQ",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for K = Int64(v_["1"]):Int64(v_["LGTR"])
            ig = ig_["R"*string(K)]
            arrset(pbm.grftype,ig,"gSQ")
            for I = Int64(v_["1"]):Int64(v_["LGSG"])
                ig = ig_["R"*string(K)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["PROD"*string(K)*","*string(I)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SXR2-MN-61-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "ePR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        NEGIDX = pbm.elpar[iel_][1]<=0.0
        if NEGIDX
            SCAL = 0.0
        end
        if !NEGIDX
            SCAL = 1.0
        end
        f_   = SCAL*EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SCAL*EV_[2]
            g_[2] = SCAL*EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = SCAL
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

    elseif action == "gSQ"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^2
        if nargout>1
            g_ = 2*GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2
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

