function VANDERM3(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A nonlinear equation problem, subject to monotonicity constraints.
#    The Jacobian is a dense Vandermonde matrix.
# 
#    Problems VANDERM1, VANDERM2, VANDERM3 and VANDERM4 differ by the rhs
#    of the equation.  They are increasingly degenerate.
# 
#    The problem is non-convex.
# 
#    Source:
#    A. Neumaier, private communication, 1991.
# 
#    SIF input: Ph. L. Toint, May 1993.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-NOR2-AN-V-V"
# 
#    Size of the system (N must be even)
# 
#       Alternative values for the SIF file parameters:
# IE N                   2              $-PARAMETER
# IE N                   3              $-PARAMETER
# IE N                   4              $-PARAMETER
# IE N                   5              $-PARAMETER
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "VANDERM3"

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
        v_["1"] = 1
        v_["2"] = 2
        v_["RN"] = Float64(v_["N"])
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["TMP"] = v_["RI"]/v_["RN"]
            v_["AL"*string(Int64(v_["I-1"]))] = v_["TMP"]
            v_["AL"*string(I)] = v_["TMP"]
        end
        v_["A"*string(Int64(v_["N"]))] = 1.0
        v_["A"*string(Int64(v_["1"]))] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["A"*string(Int64(v_["1"]))]  = (
                  v_["AL"*string(I)]+v_["A"*string(Int64(v_["1"]))])
        end
        for K = Int64(v_["2"]):Int64(v_["N"])
            v_["RK"] = Float64(K)
            v_["A"*string(K)] = 0.0
            for I = Int64(v_["1"]):Int64(v_["N"])
                v_["LOGAL"] = log(v_["AL"*string(I)])
                v_["KLOGAL"] = v_["LOGAL"]*v_["RK"]
                v_["ALK"] = exp(v_["KLOGAL"])
                v_["A"*string(K)] = v_["A"*string(K)]+v_["ALK"]
            end
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("E"*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(Int64(v_["1"])))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for K = Int64(v_["2"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("E"*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(K))
        end
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("M"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"M"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
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
        for K = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["E"*string(K)]] = Float64(v_["A"*string(K)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TMP"] = v_["RI-1"]/v_["RN"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["TMP"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePOWER", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PWR")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["2"]):Int64(v_["N"])
            v_["RK"] = Float64(K)
            for I = Int64(v_["1"]):Int64(v_["N"])
                ename = "E"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePOWER")
                arrset(ielftype,ie,iet_["ePOWER"])
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="PWR",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["RK"]))
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
        ig = ig_["E"*string(Int64(v_["1"]))]
        arrset(pbm.grftype,ig,"gL2")
        for K = Int64(v_["2"]):Int64(v_["N"])
            ig = ig_["E"*string(K)]
            arrset(pbm.grftype,ig,"gL2")
            for I = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["E"*string(K)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-AN-V-V"
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

    elseif action == "ePOWER"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]^pbm.elpar[iel_][1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[1]^(pbm.elpar[iel_][1]-1.0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1]  = (
                      pbm.elpar[iel_][1]*(pbm.elpar[iel_][1]-1.0)*EV_[1]^(pbm.elpar[iel_][1]-2.0))
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

