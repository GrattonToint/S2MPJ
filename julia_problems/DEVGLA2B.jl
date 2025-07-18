function DEVGLA2B(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DEVGLA2B
#    *********
# 
#    SCIPY global optimization benchmark example DeVilliersGlasser02
# 
#    Fit: y  = x_1 x_2^t tanh ( t x_3 + sin( t x_4 ) ) cos( t e^x_5 )  +  e
# 
#    version with box-constrained feasible region
# 
#    Source:  Problem from the SCIPY benchmark set
#      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
#              benchmarks/go_benchmark_functions
# 
#    SIF input: Nick Gould, Jan 2020
# 
#    classification = "C-CSBR2-MN-5-0"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DEVGLA2B"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling DEVGLA2B.")
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
        v_["M"] = 16
        v_["N"] = 5
        v_["1"] = 1
        v_["A"] = 1.27
        v_["LNA"] = log(v_["A"])
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["RIM1"] = -1.0+v_["RI"]
            v_["T"] = 0.1*v_["RIM1"]
            v_["T"*string(I)] = v_["T"]
            v_["TLNA"] = v_["T"]*v_["LNA"]
            v_["AT"] = exp(v_["TLNA"])
            v_["TP"] = 3.012*v_["T"]
            v_["TP2"] = 2.13*v_["T"]
            v_["STP2"] = sin(v_["TP2"])
            v_["TPA"] = v_["TP"]+v_["STP2"]
            v_["HTPA"] = tanh(v_["TPA"])
            v_["EC"] = exp(0.507)
            v_["ECT"] = v_["EC"]*v_["T"]
            v_["CECT"] = cos(v_["ECT"])
            v_["P"] = v_["AT"]*v_["HTPA"]
            v_["PP"] = v_["P"]*v_["CECT"]
            v_["PPP"] = 53.81*v_["PP"]
            v_["Y"*string(I)] = v_["PPP"]
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
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = fill(1.0,pb.n)
        pb.xupper = fill(60.0,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(20.0)
        pb.x0[ix_["X2"]] = Float64(2.0)
        pb.x0[ix_["X3"]] = Float64(2.0)
        pb.x0[ix_["X4"]] = Float64(2.0)
        pb.x0[ix_["X5"]] = Float64(0.2)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eDG2", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDG2")
            arrset(ielftype,ie,iet_["eDG2"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0),Float64(60.0),nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0),Float64(60.0),nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0),Float64(60.0),nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0),Float64(60.0),nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0),Float64(60.0),nothing)
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["T"*string(I)]))
        end
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSBR2-MN-5-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eDG2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        X2T = EV_[2]^pbm.elpar[iel_][1]
        F2 = X2T
        F2X2 = pbm.elpar[iel_][1]*EV_[2]^(pbm.elpar[iel_][1]-1.0e0)
        F2X2X2  = (
              pbm.elpar[iel_][1]*(pbm.elpar[iel_][1]-1.0e0)*EV_[2]^(pbm.elpar[iel_][1]-2.0e0))
        X3T = EV_[3]*pbm.elpar[iel_][1]
        X4T = EV_[4]*pbm.elpar[iel_][1]
        SINX4T = sin(X4T)
        COSX4T = cos(X4T)
        A = X3T+SINX4T
        AX3 = pbm.elpar[iel_][1]
        AX4 = pbm.elpar[iel_][1]*COSX4T
        AX4X4 = -pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*SINX4T
        F = tanh(A)
        FA = 1.0/(cosh(A))^2
        FAA = -2.0*FA*F
        F3 = F
        F3X3 = FA*AX3
        F3X4 = FA*AX4
        F3X3X3 = FAA*AX3*AX3
        F3X3X4 = FAA*AX3*AX4
        F3X4X4 = FA*AX4X4+FAA*AX4*AX4
        EX5 = exp(EV_[5])
        TEX5 = pbm.elpar[iel_][1]*EX5
        STEX5 = sin(TEX5)
        CTEX5 = cos(TEX5)
        F4 = CTEX5
        F4X5 = -STEX5*TEX5
        F4X5X5 = -STEX5*TEX5-CTEX5*TEX5*TEX5
        f_   = EV_[1]*F2*F3*F4
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = F2*F3*F4
            g_[2] = EV_[1]*F2X2*F3*F4
            g_[3] = EV_[1]*F2*F3X3*F4
            g_[4] = EV_[1]*F2*F3X4*F4
            g_[5] = EV_[1]*F2*F3*F4X5
            if nargout>2
                H_ = zeros(Float64,5,5)
                H_[1,2] = F2X2*F3*F4
                H_[2,1] = H_[1,2]
                H_[1,3] = F2*F3X3*F4
                H_[3,1] = H_[1,3]
                H_[1,4] = F2*F3X4*F4
                H_[4,1] = H_[1,4]
                H_[1,5] = F2*F3*F4X5
                H_[5,1] = H_[1,5]
                H_[2,2] = EV_[1]*F2X2X2*F3*F4
                H_[2,3] = EV_[1]*F2X2*F3X3*F4
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*F2X2*F3X4*F4
                H_[4,2] = H_[2,4]
                H_[2,5] = EV_[1]*F2X2*F3*F4X5
                H_[5,2] = H_[2,5]
                H_[3,3] = EV_[1]*F2*F3X3X3*F4
                H_[3,4] = EV_[1]*F2*F3X3X4*F4
                H_[4,3] = H_[3,4]
                H_[3,5] = EV_[1]*F2*F3X3*F4X5
                H_[5,3] = H_[3,5]
                H_[4,4] = EV_[1]*F2*F3X4X4*F4
                H_[4,5] = EV_[1]*F2*F3X4*F4X5
                H_[5,4] = H_[4,5]
                H_[5,5] = EV_[1]*F2*F3*F4X5X5
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

