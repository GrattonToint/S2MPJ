function DTOC3(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC3
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 1 control variable and 2 state variables.
# 
#    The problem is convex.
# 
#    Sources: problem 3 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    D.P. Bertsekas,
#    "Projected Newton methods for optimization problems with simple
#    constraints", 
#    SIAM J. Control and Optimization 20, pp. 221-246, 1982.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "C-CQLR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has 3N-1  variables (of which 2 are fixed),
#    and 2(N-1) constraints
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n=   29,m= 18 original value
# IE N                   50             $-PARAMETER  n=  149,m= 98
# IE N                   100            $-PARAMETER  n=  299,m=198
# IE N                   500            $-PARAMETER  n= 1499,m=998
# IE N                   1000           $-PARAMETER  n= 2999,m=1998
# IE N                   1500           $-PARAMETER  n= 4499,m=2998
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DTOC3"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling DTOC3.")
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
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   5000           $-PARAMETER  n=14999,m=9998
        v_["N-1"] = -1+v_["N"]
        v_["1"] = 1
        v_["2"] = 2
        v_["RN"] = Float64(v_["N"])
        v_["S"] = 1.0/v_["RN"]
        v_["2/S"] = 2.0/v_["S"]
        v_["-S"] = -1.0*v_["S"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("X"*string(T),ix_)
            arrset(pb.xnames,iv,"X"*string(T))
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["2"])
                iv,ix_,_ = s2mpj_ii("Y"*string(T)*","*string(I),ix_)
                arrset(pb.xnames,iv,"Y"*string(T)*","*string(I))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("O"*string(T),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["2/S"]))
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            v_["T+1"] = 1+T
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(Int64(v_["T+1"]))*","*string(Int64(v_["1"]))])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(T)*","*string(Int64(v_["1"]))])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["1"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["1"])))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(T)*","*string(Int64(v_["2"]))])
            push!(valA,Float64(v_["S"]))
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["2"])))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(Int64(v_["T+1"]))*","*string(Int64(v_["2"]))])
            push!(valA,Float64(-1.0))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(T)*","*string(Int64(v_["2"]))])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["2"])))
            push!(irA,ig)
            push!(icA,ix_["Y"*string(T)*","*string(Int64(v_["1"]))])
            push!(valA,Float64(v_["-S"]))
            ig,ig_,_ = s2mpj_ii("TT"*string(T)*","*string(Int64(v_["2"])),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T)*","*string(Int64(v_["2"])))
            push!(irA,ig)
            push!(icA,ix_["X"*string(T)])
            push!(valA,Float64(v_["S"]))
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 15.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]] = 15.0
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]] = 5.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]] = 5.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(15.0))
        pb.x0[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(5.0))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"YY")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["2"]):Int64(v_["N"])
            ename = "Y1SQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "Y"*string(T)*","*string(Int64(v_["1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "Y2SQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "Y"*string(T)*","*string(Int64(v_["2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "XSQ"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            v_["T+1"] = 1+T
            ig = ig_["O"*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Y1SQ"*string(Int64(v_["T+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(2.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["Y2SQ"*string(Int64(v_["T+1"]))])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(6.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  10)      224.590381002
# LO SOLUTION(  50)      233.278523083
# LO SOLUTION( 100)      234.286202920
# LO SOLUTION( 500)      235.084407947
# LO SOLUTION(1000)      235.182824435
# LO SOLUTION(5000)      235.154640099
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CQLR2-AN-V-V"
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

