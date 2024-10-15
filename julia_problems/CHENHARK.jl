function CHENHARK(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHENHARK
#    --------
# 
#    A bound-constrained version the Linear Complementarity problem
# 
#    Find x such that w = M x + q, x and w nonnegative and x^T w = 0,
#    where
# 
#    M = (  6   -4   1   0  ........ 0 ) 
#        ( -4    6  -4   1  ........ 0 )
#        (  1   -4   6  -4  ........ 0 )
#        (  0    1  -4   6  ........ 0 )  
#           ..........................
#        (  0   ........... 0  1 -4  6 )
# 
#    and q is given.
# 
#    Source: 
#    B. Chen and P. T. Harker,
#    SIMAX 14 (1993) 1168-1190
# 
#    SDIF input: Nick Gould, November 1993.
# 
#    classification = "C-QBR2-AN-V-V"
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER     original value
# IE N                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CHENHARK"

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
# IE N                   10000          $-PARAMETER
# IE N                   50000          $-PARAMETER
# IE NFREE               5              $-PARAMETER
# IE NFREE               50             $-PARAMETER
# IE NFREE               500            $-PARAMETER     original value
# IE NFREE               2500           $-PARAMETER
        if nargin<2
            v_["NFREE"] = Int64(5);  #  SIF file default value
        else
            v_["NFREE"] = Int64(args[2]);
        end
# IE NFREE               5000           $-PARAMETER
# IE NFREE               10000          $-PARAMETER
# IE NDEGEN              2              $-PARAMETER
# IE NDEGEN              20             $-PARAMETER
# IE NDEGEN              200            $-PARAMETER     original value
# IE NDEGEN              500            $-PARAMETER
        if nargin<3
            v_["NDEGEN"] = Int64(2);  #  SIF file default value
        else
            v_["NDEGEN"] = Int64(args[3]);
        end
# IE NDEGEN              1000           $-PARAMETER
# IE NDEGEN              2000           $-PARAMETER
        v_["-1"] = -1
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N-1"] = -1+v_["N"]
        v_["N+1"] = 1+v_["N"]
        v_["N+2"] = 2+v_["N"]
        v_["NFREE+1"] = 1+v_["NFREE"]
        v_["NF+ND"] = v_["NFREE"]+v_["NDEGEN"]
        v_["NF+ND+1"] = 1+v_["NF+ND"]
        v_["X"*string(Int64(v_["-1"]))] = 0.0
        v_["X"*string(Int64(v_["0"]))] = 0.0
        for I = Int64(v_["1"]):Int64(v_["NFREE"])
            v_["X"*string(I)] = 1.0
        end
        for I = Int64(v_["NFREE+1"]):Int64(v_["N+2"])
            v_["X"*string(I)] = 0.0
        end
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
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("Q"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-2.0)
        end
        ig,ig_,_ = s2mpj_ii("Q"*string(Int64(v_["0"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("Q"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("Q"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X"*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("Q"*string(Int64(v_["N+1"])),ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["NF+ND"])
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            v_["I-1"] = -1+I
            v_["I-2"] = -2+I
            v_["Q1"] = -6.0*v_["X"*string(I)]
            v_["Q2"] = 4.0*v_["X"*string(Int64(v_["I+1"]))]
            v_["Q3"] = 4.0*v_["X"*string(Int64(v_["I-1"]))]
            v_["Q4"] = -1.0*v_["X"*string(Int64(v_["I+2"]))]
            v_["Q5"] = -1.0*v_["X"*string(Int64(v_["I-2"]))]
            v_["Q"] = v_["Q1"]+v_["Q2"]
            v_["Q"] = v_["Q"]+v_["Q3"]
            v_["Q"] = v_["Q"]+v_["Q4"]
            v_["Q"] = v_["Q"]+v_["Q5"]
            ig,ig_,_ = s2mpj_ii("L",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["Q"])
        end
        for I = Int64(v_["NF+ND+1"]):Int64(v_["N"])
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            v_["I-1"] = -1+I
            v_["I-2"] = -2+I
            v_["Q1"] = -6.0*v_["X"*string(I)]
            v_["Q2"] = 4.0*v_["X"*string(Int64(v_["I+1"]))]
            v_["Q3"] = 4.0*v_["X"*string(Int64(v_["I-1"]))]
            v_["Q4"] = -1.0*v_["X"*string(Int64(v_["I+2"]))]
            v_["Q5"] = -1.0*v_["X"*string(Int64(v_["I-2"]))]
            v_["Q"] = v_["Q1"]+v_["Q2"]
            v_["Q"] = v_["Q"]+v_["Q3"]
            v_["Q"] = v_["Q"]+v_["Q4"]
            v_["Q"] = v_["Q"]+v_["Q5"]
            v_["Q"] = 1.0+v_["Q"]
            ig,ig_,_ = s2mpj_ii("L",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["Q"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.5)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gHALFL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N+1"])
            ig = ig_["Q"*string(I)]
            arrset(pbm.grftype,ig,"gHALFL2")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 1.0
#    Solution
# LO SOLTN               -0.5
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-QBR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gHALFL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= 5.0e-1*GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 1.0e+0
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

