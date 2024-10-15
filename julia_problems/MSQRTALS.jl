function MSQRTALS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MSQRTALS
#    *********
# 
#    The dense matrix square root problem by Nocedal and Liu (Case 0).
# 
#    This is a least-squares variant of problem MSQRTA.
# 
#    Source:  problem 201 (p. 93) in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-SUR2-AN-V-V"
# 
#    Dimension of the matrix
# 
#       Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER n = 4     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   10             $-PARAMETER n = 100
# IE P                   23             $-PARAMETER n = 529
# IE P                   32             $-PARAMETER n = 1024
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MSQRTALS"

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
            v_["P"] = Int64(5);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
# IE P                   70             $-PARAMETER n = 4900
        v_["N"] = v_["P"]*v_["P"]
        v_["1"] = 1
        v_["K"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                v_["K"] = 1.0+v_["K"]
                v_["K2"] = v_["K"]*v_["K"]
                v_["B"*string(I)*","*string(J)] = sin(v_["K2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                v_["A"*string(I)*","*string(J)] = 0.0
                for K = Int64(v_["1"]):Int64(v_["P"])
                    v_["PROD"]  = (
                          v_["B"*string(I)*","*string(K)]*v_["B"*string(K)*","*string(J)])
                    v_["A"*string(I)*","*string(J)] = (v_["A"*string(I)*","*string(J)]+
                         v_["PROD"])
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                ig,ig_,_ = s2mpj_ii("G"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                pbm.gconst[ig_["G"*string(I)*","*string(J)]]  = (
                      Float64(v_["A"*string(I)*","*string(J)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        v_["K"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                v_["K"] = 1.0+v_["K"]
                v_["K2"] = v_["K"]*v_["K"]
                v_["SK2"] = sin(v_["K2"])
                v_["-4SK2/5"] = -0.8*v_["SK2"]
                v_["XIJ"] = v_["B"*string(I)*","*string(J)]+v_["-4SK2/5"]
                pb.x0[ix_["X"*string(I)*","*string(J)]] = Float64(v_["XIJ"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"XIT")
        loaset(elftv,it,2,"XTJ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                for T = Int64(v_["1"]):Int64(v_["P"])
                    ename = "E"*string(I)*","*string(J)*","*string(T)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                    vname = "X"*string(I)*","*string(T)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="XIT",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(T)*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="XTJ",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
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
        for I = Int64(v_["1"]):Int64(v_["P"])
            for J = Int64(v_["1"]):Int64(v_["P"])
                for T = Int64(v_["1"]):Int64(v_["P"])
                    ig = ig_["G"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)*","*string(T)])
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SUR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
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

