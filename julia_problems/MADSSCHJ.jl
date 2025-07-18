function MADSSCHJ(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MADSSCHJ
#    *********
# 
#    A nonlinear minmax problem with variable dimension.
#    The Jacobian of the constraints is dense.
# 
#    Source:
#    K. Madsen and H. Schjaer-Jacobsen,
#    "Linearly Constrained Minmax Optimization",
#    Mathematical Programming 14, pp. 208-223, 1978.
# 
#    SIF input: Ph. Toint, August 1993.
# 
#    classification = "C-CLQR2-AN-V-V"
# 
#    N is the number of variables - 1, and must be even and at least 4.
#    The number of inequality constraints is 2*N - 2.
# 
#       Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER  n=  5, m=  6
# IE N                   10             $-PARAMETER  n= 11, m= 18  original value
# IE N                   20             $-PARAMETER  n= 21, m= 38
# IE N                   30             $-PARAMETER  n= 31, m= 58
# IE N                   40             $-PARAMETER  n= 41, m= 78
# IE N                   50             $-PARAMETER  n= 51, m= 98
# IE N                   60             $-PARAMETER  n= 61, m=118
# IE N                   70             $-PARAMETER  n= 71, m=138
# IE N                   80             $-PARAMETER  n= 81, m=158
# IE N                   90             $-PARAMETER  n= 91, m=178
# IE N                   100            $-PARAMETER  n=101, m=198
# IE N                   200            $-PARAMETER  n=201, m=398
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MADSSCHJ"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling MADSSCHJ.")
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
            v_["N"] = Int64(4);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["N-1"] = -1+v_["N"]
        v_["2N"] = v_["N"]+v_["N"]
        v_["M"] = -2+v_["2N"]
        v_["M-1"] = -1+v_["M"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("Z",ix_)
        arrset(pb.xnames,iv,"Z")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["Z"])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C1")
        push!(irA,ig)
        push!(icA,ix_["Z"])
        push!(valA,Float64(1.0))
        for I = Int64(v_["2"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("C1",ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C1")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
        end
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C2")
        push!(irA,ig)
        push!(icA,ix_["Z"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["3"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("C2",ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C2")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
        end
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C3")
        push!(irA,ig)
        push!(icA,ix_["Z"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["3"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("C3",ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C3")
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
        end
        for K = Int64(v_["4"]):Int64(v_["2"]):Int64(v_["M-1"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            v_["J"] = trunc(Int,(v_["K+2"]/v_["2"]))
            v_["J-1"] = -1+v_["J"]
            v_["J+1"] = 1+v_["J"]
            ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(K))
            push!(irA,ig)
            push!(icA,ix_["Z"])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["K+1"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["K+1"])))
            push!(irA,ig)
            push!(icA,ix_["Z"])
            push!(valA,Float64(1.0))
            for I = Int64(v_["1"]):Int64(v_["J-1"])
                ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(K))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)])
                push!(valA,Float64(-1.0))
                ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["K+1"])),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(Int64(v_["K+1"])))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)])
                push!(valA,Float64(-1.0))
            end
            for I = Int64(v_["J+1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(K))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)])
                push!(valA,Float64(-1.0))
                ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["K+1"])),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(Int64(v_["K+1"])))
                push!(irA,ig)
                push!(icA,ix_["X"*string(I)])
                push!(valA,Float64(-1.0))
            end
        end
        ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["M"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["M"])))
        push!(irA,ig)
        push!(icA,ix_["Z"])
        push!(valA,Float64(1.0))
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["M"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["M"])))
            push!(irA,ig)
            push!(icA,ix_["X"*string(I)])
            push!(valA,Float64(-1.0))
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
        for K = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["C"*string(K)]] = Float64(-1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(10.0),pb.n)
        pb.x0[ix_["Z"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "XSQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(10.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["C"*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C"*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C"*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2.0))
        for K = Int64(v_["4"]):Int64(v_["2"]):Int64(v_["M-1"])
            v_["K+1"] = 1+K
            v_["K+2"] = 2+K
            v_["J"] = trunc(Int,(v_["K+2"]/v_["2"]))
            ig = ig_["C"*string(K)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["J"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["C"*string(Int64(v_["K+1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["J"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-2.0))
        end
        ig = ig_["C"*string(Int64(v_["M"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(4)             -2.6121094144
# LO SOLTN(10)            -12.814452425
# LO SOLTN(20)            -49.869888156
# LO SOLTN(30)            -111.93545559
# LO SOLTN(40)            -199.00371592
# LO SOLTN(50)            -311.07308068
# LO SOLTN(60)            -448.14300524
# LO SOLTN(70)            -610.21325256
# LO SOLTN(80)            -797.28370289
# LO SOLTN(90)            -1009.3542892
# LO SOLTN(100)           -1246.4249710
# LO SOLTN(200)           -4992.1339031
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLQR2-AN-V-V"
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

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]+EV_[1]
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

