function DEVGLA1NE(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEVGLA1NE
#    *********
# 
#    SCIPY global optimization benchmark example DeVilliersGlasser01
# 
#    Fit: y  = x_1 x_2^t sin( t x_3 + x_4 )  +  e
# 
#    Source:  Problem from the SCIPY benchmark set
#      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
#              benchmarks/go_benchmark_functions
# 
#    Nonlinear-equation formulation of DEVGLA1.SIF
# 
#    SIF input: Nick Gould, Jan 2020
# 
#    classification = "C-CNOR2-MN-4-24"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DEVGLA1NE"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling DEVGLA1NE.")
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
        v_["M"] = 24
        v_["N"] = 4
        v_["1"] = 1
        v_["A"] = 1.371
        v_["LNA"] = log(v_["A"])
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["RIM1"] = -1.0+v_["RI"]
            v_["T"] = 0.1*v_["RIM1"]
            v_["T"*string(I)] = v_["T"]
            v_["TLNA"] = v_["T"]*v_["LNA"]
            v_["AT"] = exp(v_["TLNA"])
            v_["TP"] = 3.112*v_["T"]
            v_["TPA"] = 1.761+v_["TP"]
            v_["STPA"] = sin(v_["TPA"])
            v_["P"] = v_["AT"]*v_["STPA"]
            v_["PP"] = 60.137*v_["P"]
            v_["Y"*string(I)] = v_["PP"]
        end
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(2.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eDG1", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDG1")
            arrset(ielftype,ie,iet_["eDG1"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-MN-4-24"
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

    elseif action == "eDG1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = pbm.elpar[iel_][1]*EV_[3]+EV_[4]
        SINA = sin(A)
        COSA = cos(A)
        X2T = EV_[2]^pbm.elpar[iel_][1]
        DX2T = pbm.elpar[iel_][1]*EV_[2]^(pbm.elpar[iel_][1]-1.0e0)
        D2X2T  = (
              pbm.elpar[iel_][1]*(pbm.elpar[iel_][1]-1.0e0)*EV_[2]^(pbm.elpar[iel_][1]-2.0e0))
        X1X2T = EV_[1]*X2T
        f_   = X1X2T*SINA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = X2T*SINA
            g_[2] = EV_[1]*DX2T*SINA
            g_[3] = X1X2T*COSA*pbm.elpar[iel_][1]
            g_[4] = X1X2T*COSA
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = DX2T*SINA
                H_[2,1] = H_[1,2]
                H_[1,3] = X2T*COSA*pbm.elpar[iel_][1]
                H_[3,1] = H_[1,3]
                H_[1,4] = X2T*COSA
                H_[4,1] = H_[1,4]
                H_[2,2] = EV_[1]*D2X2T*SINA
                H_[2,3] = EV_[1]*DX2T*COSA*pbm.elpar[iel_][1]
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*DX2T*COSA
                H_[4,2] = H_[2,4]
                H_[3,3] = -X1X2T*SINA*pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
                H_[3,4] = -X1X2T*SINA*pbm.elpar[iel_][1]
                H_[4,3] = H_[3,4]
                H_[4,4] = -X1X2T*SINA
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

