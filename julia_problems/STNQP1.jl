function STNQP1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : STNQP1
#    *********
# 
#    A non-convex quadratic program with some structure.
# 
#    The objective function is of the form
#       sum (i=0,n) x_i^2 - 0.5 sum (l=1,n/p) sum(i=1,p) sum(k;i) x_{k+l}^2,
#    where n = 2^p and (k;i) means k takes the values of the first i powers of 2
#    eg, (k:3) = {k = {1,2,4}} and (k:7) = {k = {1,2,4,8,16,32}}.
#    There are equality constraints of the form
#    
#       sum(k;i) x_{k+l-1} = i, where l=1,n/p,2 and i=1,p.
#    Finally, there are simple bounds
#          2 <= x_i, y_i <= 2    (i=0,n).
# 
#    SIF input: Nick Gould, May 1996
# 
#    classification = "C-QLR2-AN-V-V"
# 
#    There will be 2**p + 1 variables
# 
#       Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER n = 5
# IE P                   4              $-PARAMETER n = 17
# IE P                   6              $-PARAMETER n = 65
# IE P                   8              $-PARAMETER n = 257
# IE P                   10             $-PARAMETER n = 1025
# IE P                   12             $-PARAMETER n = 4097     original value
# IE P                   13             $-PARAMETER n = 8193
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "STNQP1"

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
# IE P                   14             $-PARAMETER n = 16395
# IE P                   15             $-PARAMETER n = 32769
# IE P                   16             $-PARAMETER n = 65537
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N"] = 1
        for I = Int64(v_["1"]):Int64(v_["P"])
            v_["N"] = v_["N"]*v_["2"]
        end
        v_["N/P"] = trunc(Int,(v_["N"]/v_["P"]))
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for L = Int64(v_["1"]):Int64(v_["N/P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                v_["K"] = v_["1"]
                for J = Int64(v_["1"]):Int64(I)
                    v_["K+L"] = v_["K"]+L
                    ig,ig_,_ = s2mpj_ii("N"*string(I)*","*string(L),ig_)
                    arrset(gtype,ig,"<>")
                    iv = ix_["X"*string(Int64(v_["K+L"]))]
                    pbm.A[ig,iv] += Float64(1.0)
                    v_["K"] = v_["K"]*v_["2"]
                end
            end
        end
        for L = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["N/P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                v_["K"] = v_["1"]
                for J = Int64(v_["1"]):Int64(I)
                    v_["K-1"] = v_["K"]-v_["1"]
                    v_["K+L-1"] = v_["K-1"]+L
                    ig,ig_,_ = s2mpj_ii("E"*string(I)*","*string(L),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"E"*string(I)*","*string(L))
                    iv = ix_["X"*string(Int64(v_["K+L-1"]))]
                    pbm.A[ig,iv] += Float64(1.0)
                    v_["K"] = v_["K"]*v_["2"]
                end
            end
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
        for L = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["N/P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                v_["RI"] = Float64(I)
                pbm.gconst[ig_["E"*string(I)*","*string(L)]] = Float64(v_["RI"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["0"]):Int64(v_["N"])
            pb.xlower[ix_["X"*string(I)]] = -2.0
            pb.xupper[ix_["X"*string(I)]] = 2.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gPSQR",igt_)
        it,igt_,_ = s2mpj_ii("gPSQR",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            ig = ig_["O"*string(I)]
            arrset(pbm.grftype,ig,"gPSQR")
            posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
            loaset(pbm.grpar,ig,posgp,Float64(1.0))
        end
        for L = Int64(v_["1"]):Int64(v_["N/P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                ig = ig_["N"*string(I)*","*string(L)]
                arrset(pbm.grftype,ig,"gPSQR")
                posgp = findfirst(x->x=="P",grftp[igt_[pbm.grftype[ig]]])
                loaset(pbm.grpar,ig,posgp,Float64(-0.5))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION            -1.361565E+5   $ (P=12)
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

    elseif action == "gPSQR"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= pbm.grpar[igr_][1]*GVAR_*GVAR_
        if nargout>1
            g_ = 2.0*pbm.grpar[igr_][1]*GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0*pbm.grpar[igr_][1]
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

