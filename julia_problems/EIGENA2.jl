function EIGENA2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : EIGENA2
#    --------
# 
#    Solving symmetric eigenvalue problems as systems of
#    nonlinear equations.
# 
#    The problem is, given a symmetric matrix A, to find an orthogonal
#    matrix Q and diagonal matrix D such that A Q(T) = Q(T) D.
# 
#    Example A: a diagonal matrix with eigenvales 1, .... , N.
# 
#    Source:  An idea by Nick Gould
# 
#                Constrained optimization version 2.
# 
#    SIF input: Nick Gould, Nov 1992.
# 
#    classification = "C-QQR2-AN-V-V"
# 
#    The dimension of the matrix.
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "EIGENA2"

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
            v_["N"] = Int64(2);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   100            $-PARAMETER
        v_["1"] = 1
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                v_["A"*string(I)*","*string(J)] = 0.0
            end
            v_["RJ"] = Float64(J)
            v_["A"*string(J)*","*string(J)] = v_["RJ"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("D"*string(J),ix_)
            arrset(pb.xnames,iv,"D"*string(J))
            for I = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("Q"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Q"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(J)
                ig,ig_,_ = s2mpj_ii("O"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"O"*string(I)*","*string(J))
            end
            for I = Int64(v_["1"]):Int64(v_["N"])
                for K = Int64(v_["1"]):Int64(v_["N"])
                    v_["-AIK"] = -1.0*v_["A"*string(I)*","*string(K)]
                    ig,ig_,_ = s2mpj_ii("E"*string(I)*","*string(J),ig_)
                    arrset(gtype,ig,"<>")
                    iv = ix_["Q"*string(J)*","*string(K)]
                    pbm.A[ig,iv] += Float64(v_["-AIK"])
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
        for J = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["O"*string(J)*","*string(J)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["RJ"] = Float64(J)
            pb.x0[ix_["D"*string(J)]] = Float64(1.0)
            pb.x0[ix_["Q"*string(J)*","*string(J)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PROD", iet_)
        loaset(elftv,it,1,"Q1")
        loaset(elftv,it,2,"Q2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ename = "E"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"en2PROD")
                    arrset(ielftype,ie,iet_["en2PROD"])
                end
                vname = "Q"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="Q1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "D"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="Q2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for I = Int64(v_["1"]):Int64(J)
                for K = Int64(v_["1"]):Int64(v_["N"])
                    ename = "O"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,newelt = s2mpj_ii(ename,ie_)
                    if newelt > 0
                        arrset(pbm.elftype,ie,"en2PROD")
                        arrset(ielftype,ie,iet_["en2PROD"])
                    end
                    vname = "Q"*string(K)*","*string(I)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                    posev = findfirst(x->x=="Q1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "Q"*string(K)*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                    posev = findfirst(x->x=="Q2",elftv[ielftype[ie]])
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
        for J = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["E"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
            for I = Int64(v_["1"]):Int64(J)
                for K = Int64(v_["1"]):Int64(v_["N"])
                    ig = ig_["O"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["O"*string(I)*","*string(J)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
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
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QQR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en2PROD"

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
                H_[1,2] = 1.0e+0
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
                H_ = 2.0e+0
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

