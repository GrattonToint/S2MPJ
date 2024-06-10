function DIAGPQT(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DIAGPQT
#    *********
#    A variable dimension convex quadratic problem
#    with eigenvalues clusetered towards the top of the spectrum
# 
#    lambda_i = - i^2/n + n + 1/n, i = 1, ... , n
# 
#    Source: simple test for GALAHAD gltr/glrt
# 
#    SIF input: Nick Gould, Feb 2019
# 
#    classification = "QBR2-AN-V-0"
# 
#    Number of variables (variable)
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DIAGPQT"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER     original value
# IE N                   1000           $-PARAMETER
# IE N                   5000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# IE N                   100000         $-PARAMETER
# IE N                   1000000        $-PARAMETER
        v_["1"] = 1
        v_["ONE"] = 1.0
        v_["RN"] = Float64(v_["N"])
        v_["1/RN"] = v_["ONE"]/v_["RN"]
        v_["SHIFT"] = v_["RN"]+v_["1/RN"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-100000.0,pb.n)
        pb.xupper = fill(1000000.0,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["RI2"] = v_["RI"]*v_["RI"]
            v_["RI2/RN"] = v_["RI2"]/v_["RN"]
            v_["-RI2/RN"] = -1.0*v_["RI2/RN"]
            v_["H"] = v_["-RI2/RN"]+v_["SHIFT"]
            ix1 = ix_["X"*string(I)]
            ix2 = ix_["X"*string(I)]
            pbm.H[ix1,ix2] = Float64(v_["H"])+pbm.H[ix1,ix2]
            pbm.H[ix2,ix1] = pbm.H[ix1,ix2]
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        Hsave = pbm.H[ 1:pb.n, 1:pb.n ]
        pbm.H = Hsave
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "QBR2-AN-V-0"
        return pb, pbm

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

