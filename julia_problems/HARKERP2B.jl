function HARKERP2B(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HARKERP2B
#    --------
# 
#    A bound-constrained version of a Linear Complementarity problem
#    posed by Harker and Pang.
# 
#    A test version without group scaling.
# 
#    Source: 
#    P. T. Harker and J.-S. Pang,
#    "A damped Newton method for the linear complementarity problem",
#    in 'Allgower and Georg: Computational solution of nonlinear
#    systems of equations', AMS lectures in Applied Mathematics 26,
#    AMS, Providence, Rhode Island, USA, pp 265-284.
# 
#    SIF input: Nick Gould, July 1993.
# 
#    classification = "QBR2-AN-V-V"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HARKERP2B"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "HARKERP2B"
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
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2x_ii("S"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(-1.0))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2x_ii("Q"*string(Int64(v_["0"])),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2x_ii("Q"*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        ig,ig_,_ = s2x_ii("Q"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(0.5))
        for J = Int64(v_["2"]):Int64(v_["N"])
            ig,ig_,_ = s2x_ii("Q"*string(J),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(0.25))
            for I = Int64(J):Int64(v_["N"])
                ig,ig_,_ = s2x_ii("Q"*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            pb.x0[ix_["X"*string(I)]] = Float64(v_["RI"])
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2x_ii("gHALFL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["S"*string(I)]
            arrset(pbm.grftype,ig,"gHALFL2")
            ig = ig_["Q"*string(I)]
            arrset(pbm.grftype,ig,"gHALFL2")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xlower =    fill(Inf,pb.n)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "QBR2-AN-V-V"
        return pb, pbm

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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
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

