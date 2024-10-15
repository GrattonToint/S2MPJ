function MINSURF(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MINSURF
#    *********
#    Variable dimension full rank linear problem
#    A version of the minimum surface problem
#    on the unit square with simple boundary conditions.
# 
#    SIF input: Ph. Toint, Jan 1991.
# 
#    classification = "C-OXR2-MY-64-0"
# 
#    Discretization parameter
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MINSURF"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["P"] = 7
        v_["1"] = 1
        v_["P+1"] = 1+v_["P"]
        v_["RP"] = Float64(v_["P"])
        v_["RPSQ"] = v_["RP"]*v_["RP"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for i = Int64(v_["1"]):Int64(v_["P+1"])
            for j = Int64(v_["1"]):Int64(v_["P+1"])
                iv,ix_,_ = s2mpj_ii("X"*string(i)*","*string(j),ix_)
                arrset(pb.xnames,iv,"X"*string(i)*","*string(j))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for i = Int64(v_["1"]):Int64(v_["P"])
            for j = Int64(v_["1"]):Int64(v_["P"])
                ig,ig_,_ = s2mpj_ii("S"*string(i)*","*string(j),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["RPSQ"]))
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for i = Int64(v_["1"]):Int64(v_["P"])
            for j = Int64(v_["1"]):Int64(v_["P"])
                pbm.gconst[ig_["S"*string(i)*","*string(j)]] = Float64(-1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        v_["2"] = 2
        for i = Int64(v_["2"]):Int64(v_["P"])
            for j = Int64(v_["2"]):Int64(v_["P"])
                pb.xlower[ix_["X"*string(i)*","*string(j)]] = -Inf
                pb.xupper[ix_["X"*string(i)*","*string(j)]] = +Inf
            end
        end
        for i = Int64(v_["1"]):Int64(v_["P+1"])
            pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(i)]] = 1.0
            pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(i)]] = 1.0
            pb.xlower[ix_["X"*string(Int64(v_["P+1"]))*","*string(i)]] = 1.0
            pb.xupper[ix_["X"*string(Int64(v_["P+1"]))*","*string(i)]] = 1.0
            pb.xlower[ix_["X"*string(i)*","*string(Int64(v_["1"]))]] = 1.0
            pb.xupper[ix_["X"*string(i)*","*string(Int64(v_["1"]))]] = 1.0
            pb.xlower[ix_["X"*string(i)*","*string(Int64(v_["P+1"]))]] = 1.0
            pb.xupper[ix_["X"*string(i)*","*string(Int64(v_["P+1"]))]] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["1"]):Int64(v_["P"])
            v_["i+1"] = 1+i
            for j = Int64(v_["1"]):Int64(v_["P"])
                v_["j+1"] = 1+j
                ename = "A"*string(i)*","*string(j)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(i)*","*string(j)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["i+1"]))*","*string(Int64(v_["j+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="W",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "B"*string(i)*","*string(j)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(i)*","*string(Int64(v_["j+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["i+1"]))*","*string(j)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="W",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQROOT",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        v_["WEIGHT"] = 0.5*v_["RPSQ"]
        for i = Int64(v_["1"]):Int64(v_["P"])
            for j = Int64(v_["1"]):Int64(v_["P"])
                ig = ig_["S"*string(i)*","*string(j)]
                arrset(pbm.grftype,ig,"gSQROOT")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(i)*","*string(j)])
                loaset(pbm.grelw,ig,posel,Float64(v_["WEIGHT"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(i)*","*string(j)])
                loaset(pbm.grelw,ig,posel,Float64(v_["WEIGHT"]))
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OXR2-MY-64-0"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
                H_ = U_'*H_*U_
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

    elseif action == "gSQROOT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SQRAL = sqrt(GVAR_)
        f_= SQRAL
        if nargout>1
            g_ = 0.5/SQRAL
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -0.25/(GVAR_*SQRAL)
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

