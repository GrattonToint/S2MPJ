function NCVXQP6(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NCVXQP6
#    *********
# 
#    A non-convex quadratic program.
# 
#    SIF input: Nick Gould, April 1995
# 
#    classification = "C-QLR2-AN-V-V"
# 
#    The number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER    original value
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "NCVXQP6"

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
            v_["N"] = Int64(50);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   100000         $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["M"] = trunc(Int,(v_["N"]/v_["4"]))
        v_["NPLUS"] = trunc(Int,(v_["N"]/v_["4"]))
        v_["NPLUS"] = v_["NPLUS"]*v_["3"]
        v_["NPLUS+1"] = 1+v_["NPLUS"]
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            v_["J"] = 2*I
            v_["J"] = -1+v_["J"]
            v_["K"] = trunc(Int,(v_["J"]/v_["N"]))
            v_["K"] = v_["K"]*v_["N"]
            v_["J"] = v_["J"]-v_["K"]
            v_["J"] = 1+v_["J"]
            iv = ix_["X"*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
            v_["J"] = 3*I
            v_["J"] = -1+v_["J"]
            v_["K"] = trunc(Int,(v_["J"]/v_["N"]))
            v_["K"] = v_["K"]*v_["N"]
            v_["J"] = v_["J"]-v_["K"]
            v_["J"] = 1+v_["J"]
            iv = ix_["X"*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("CON"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CON"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            v_["J"] = 4*I
            v_["J"] = -1+v_["J"]
            v_["K"] = trunc(Int,(v_["J"]/v_["N"]))
            v_["K"] = v_["K"]*v_["N"]
            v_["J"] = v_["J"]-v_["K"]
            v_["J"] = 1+v_["J"]
            iv = ix_["X"*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(2.0)
            v_["J"] = 5*I
            v_["J"] = -1+v_["J"]
            v_["K"] = trunc(Int,(v_["J"]/v_["N"]))
            v_["K"] = v_["K"]*v_["N"]
            v_["J"] = v_["J"]-v_["K"]
            v_["J"] = 1+v_["J"]
            iv = ix_["X"*string(Int64(v_["J"]))]
            pbm.A[ig,iv] += Float64(3.0)
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
            pbm.gconst[ig_["CON"*string(I)]] = Float64(6.0)
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["X"*string(I)]] = 0.1
            pb.xupper[ix_["X"*string(I)]] = 10.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQR",igt_)
        it,igt_,_ = s2mpj_ii("gSQR",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NPLUS"])
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gSQR")
            v_["RI"] = Float64(I)
            posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
            loaset(pbm.grpar,ig,posgp,Float64(v_["RI"]))
        end
        for I = Int64(v_["NPLUS+1"]):Int64(v_["N"])
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gSQR")
            v_["RI"] = Float64(I)
            v_["RI"] = -1.0*v_["RI"]
            posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
            loaset(pbm.grpar,ig,posgp,Float64(v_["RI"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.54861D+07   $ (n=1000)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = collect(1:length(pbm.congrps))
        pb.pbclass = "C-QLR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQR"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= 0.5*pbm.grpar[igr_][1]*GVAR_*GVAR_
        if nargout>1
            g_ = pbm.grpar[igr_][1]*GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = pbm.grpar[igr_][1]
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

