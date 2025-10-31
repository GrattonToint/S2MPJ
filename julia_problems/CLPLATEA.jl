function CLPLATEA(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CLPLATEA
#    *********
# 
#    The clamped plate problem (Strang, Nocedal, Dax)
#    The problem comes from the discretization the following problem
#    in mechanics:  a plate is clamped on one edge and loaded on the
#    opposite side.  The plate is the unit square.
# 
#    In this version of the problem, the weight WGHT is entirely put on the
#    upper right corner of the plate.
# 
#    The plate is clamped on its lower edge, by fixing the
#    corresponding variables to zero.
# 
#    Source:
#    J. Nocedal,
#    "Solving large nonlinear systems of equations arising in mechanics",
#    Proceedings of the Cocoyoc Numerical Analysis Conference, Mexico,
#    pp. 132-141, 1981.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-COXR2-MN-V-0"
# 
#    P is the number of points in one side of the unit square
#    The number of variables is P*P, of which (P-1)*(P-1) are free.
# 
#       Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 16
# IE P                   7              $-PARAMETER n = 49    original value
# IE P                   10             $-PARAMETER n = 100
# IE P                   23             $-PARAMETER n = 529
# IE P                   32             $-PARAMETER n = 1024
# IE P                   71             $-PARAMETER n = 5041
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 31 X 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CLPLATEA"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling CLPLATEA.")
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
        if nargin<1
            v_["P"] = Int64(4);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
        v_["WGHT"] = -0.1
        v_["1"] = 1
        v_["2"] = 2
        v_["P2"] = v_["P"]*v_["P"]
        v_["RP2"] = Float64(v_["P2"])
        v_["HP2"] = 0.5*v_["RP2"]
        v_["1/HP2"] = 1.0/v_["HP2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for J = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["2"]):Int64(v_["P"])
            v_["I-1"] = -1+I
            for J = Int64(v_["2"]):Int64(v_["P"])
                v_["J-1"] = -1+J
                ig,ig_,_ = s2mpj_ii("A"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(2.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(J)])
                push!(valA,Float64(1.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(Int64(v_["J-1"]))])
                push!(valA,Float64(-1.0))
                ig,ig_,_ = s2mpj_ii("B"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(2.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(J)])
                push!(valA,Float64(1.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(Int64(v_["I-1"]))*","*string(J)])
                push!(valA,Float64(-1.0))
                ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["1/HP2"]))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(J)])
                push!(valA,Float64(1.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(Int64(v_["J-1"]))])
                push!(valA,Float64(-1.0))
                ig,ig_,_ = s2mpj_ii("D"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["1/HP2"]))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)*","*string(J)])
                push!(valA,Float64(1.0))
                push!(irA,ig)
                push!(icA,ix_["X"*string(Int64(v_["I-1"]))*","*string(J)])
                push!(valA,Float64(-1.0))
            end
        end
        ig,ig_,_ = s2mpj_ii("W",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["P"]))*","*string(Int64(v_["P"]))])
        push!(valA,Float64(v_["WGHT"]))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["1"]):Int64(v_["P"])
            pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        it,igt_,_ = s2mpj_ii("gL4",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["P"])
            for J = Int64(v_["2"]):Int64(v_["P"])
                ig = ig_["A"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["B"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["C"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL4")
                ig = ig_["D"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL4")
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)            -6.5955D-03
# LO SOLTN(7)            -8.2007D-03
# LO SOLTN(10)           -9.1452D-03
# LO SOLTN(23)           -1.0974D-02
# LO SOLTN(32)           -1.1543D-02
# LO SOLTN(71)           -1.2592D-02
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-COXR2-MN-V-0"
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

    elseif action == "gL4"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_^4
        if nargout>1
            g_ = 4.0*GVAR_^3
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 12.0*GVAR_^2
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

