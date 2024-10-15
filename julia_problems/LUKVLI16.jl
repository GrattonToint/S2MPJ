function LUKVLI16(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLI16
#    *********
# 
#    Source: Problem 5.16, the chained modified HS51 problem, 
#    due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    Equality constraints changed to inequalities
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-OQR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#       Alternative values for the SIF file parameters:
# IE N                   97             $-PARAMETER
# IE N                   997            $-PARAMETER
# IE N                   9997           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKVLI16"

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
            v_["N"] = Int64(17);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   99997          $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["N-1"] = -1+v_["N"]
        v_["(N-1)/4"] = trunc(Int,(v_["N-1"]/v_["4"]))
        v_["NC"] = v_["3"]*v_["(N-1)/4"]
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
        for I = Int64(v_["1"]):Int64(v_["(N-1)/4"])
            v_["I-1"] = -1+I
            v_["J"] = v_["4"]*v_["I-1"]
            v_["J+1"] = 1+v_["J"]
            v_["J+2"] = 2+v_["J"]
            v_["J+3"] = 3+v_["J"]
            v_["J+4"] = 4+v_["J"]
            v_["J+5"] = 5+v_["J"]
            ig,ig_,_ = s2mpj_ii("OBJ1"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["J+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["J+2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("OBJ2"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["J+2"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["J+3"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("OBJ3"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["J+4"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("OBJ4"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["J+5"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for K = Int64(v_["1"]):Int64(v_["3"]):Int64(v_["NC"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            v_["K+3"] = 3+K
            v_["K+4"] = 4+K
            ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(K))
            iv = ix_["X"*string(Int64(v_["K+1"]))]
            pbm.A[ig,iv] += Float64(3.0)
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["K+1"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["K+1"])))
            iv = ix_["X"*string(Int64(v_["K+3"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["K+4"]))]
            pbm.A[ig,iv] += Float64(-2.0)
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["K+2"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["K+2"])))
            iv = ix_["X"*string(Int64(v_["K+4"]))]
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
        for I = Int64(v_["1"]):Int64(v_["(N-1)/4"])
            pbm.gconst[ig_["OBJ2"*string(I)]] = Float64(2.0)
            pbm.gconst[ig_["OBJ3"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["OBJ4"*string(I)]] = Float64(1.0)
        end
        for K = Int64(v_["1"]):Int64(v_["3"]):Int64(v_["NC"])
            pbm.gconst[ig_["C"*string(K)]] = Float64(4.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(2.5)
        end
        for I = Int64(v_["2"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.5)
        end
        for I = Int64(v_["3"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(2.0)
        end
        for I = Int64(v_["4"]):Int64(v_["4"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(-1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["1"]):Int64(v_["3"]):Int64(v_["NC"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            ename = "E"*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "X"*string(K)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["K+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            ename = "E"*string(Int64(v_["K+1"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["K+2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["K+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            ename = "E"*string(Int64(v_["K+2"]))
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X"*string(Int64(v_["K+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        it,igt_,_ = s2mpj_ii("gL4",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["(N-1)/4"])
            ig = ig_["OBJ1"*string(I)]
            arrset(pbm.grftype,ig,"gL4")
            ig = ig_["OBJ2"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OBJ3"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OBJ4"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
        end
        for K = Int64(v_["1"]):Int64(v_["3"]):Int64(v_["NC"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            ig = ig_["C"*string(K)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(K)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["C"*string(Int64(v_["K+1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["K+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["C"*string(Int64(v_["K+2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["K+2"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 2.95007E+02
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OQR2-AY-V-V"
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

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

