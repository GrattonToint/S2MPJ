function DTOC6(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DTOC6
#    *********
# 
#    This is a discrete time optimal control (DTOC) problem.  
#    The system has N time periods, 1 control variable and 1 state variable.
# 
#    The problem is convex.
# 
#    Sources: problem 6 in
#    T.F. Coleman and A. Liao,
#    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
#    Control Problems",
#    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
#    Cornell University, 1992.
# 
#    D.M. Murray and S.J. Yakowitz,
#    "The application of optimal contraol methodology to nonlinear programming
#    problems",
#    Mathematical Programming 21, pp. 331-347, 1981.
# 
#    SIF input: Ph. Toint, August 1993
# 
#    classification = "C-OOR2-AN-V-V"
# 
#    Problem variants: they are identified by the value of the parameter N.
# 
#    The problem has 2N-1  variables (of which 1 is fixed),
#    and N-1 constraints
# 
#       Alternative values for the SIF file parameters:
# IE N                   11             $-PARAMETER n =   21, m =  10
# IE N                   21             $-PARAMETER n =   41, m =  20
# IE N                   31             $-PARAMETER n =   61, m =  30
# IE N                   41             $-PARAMETER n =   81, m =  40
# IE N                   51             $-PARAMETER n =  101, m =  50
# IE N                   61             $-PARAMETER n =  121, m =  60
# IE N                   71             $-PARAMETER n =  141, m =  70
# IE N                   81             $-PARAMETER n =  161, m =  80
# IE N                   91             $-PARAMETER n =  181, m =  90
# IE N                   101            $-PARAMETER n =  201, m = 100
# IE N                   501            $-PARAMETER n = 1001, m = 500
# IE N                   1001           $-PARAMETER n = 2001, m =1000
# IE N                   5001           $-PARAMETER n =10001, m =5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DTOC6"

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
            v_["N"] = Int64(11);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["N-1"] = -1+v_["N"]
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("X"*string(T),ix_)
            arrset(pb.xnames,iv,"X"*string(T))
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("Y"*string(T),ix_)
            arrset(pb.xnames,iv,"Y"*string(T))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("OY"*string(T),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["Y"*string(T)]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(2.0))
            ig,ig_,_ = s2mpj_ii("OX"*string(T),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(T)]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(2.0))
        end
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            v_["T+1"] = 1+T
            ig,ig_,_ = s2mpj_ii("TT"*string(T),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"TT"*string(T))
            iv = ix_["Y"*string(Int64(v_["T+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["Y"*string(T)]
            pbm.A[ig,iv] += Float64(1.0)
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
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["Y"*string(Int64(v_["1"]))]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXP", iet_)
        loaset(elftv,it,1,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "E"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
            vname = "X"*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["OX"*string(T)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OY"*string(T)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["TT"*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION(  11)      19.80411526774
# LO SOLUTION(  21)      62.49481326823
# LO SOLUTION(  31)      119.0328455446
# LO SOLUTION(  41)      185.8961987565
# LO SOLUTION(  51)      261.1131573312
# LO SOLUTION(  61)      343.3995402405
# LO SOLUTION(  71)      431.8430677467
# LO SOLUTION(  81)      525.7575776948
# LO SOLUTION(  91)      624.6051839906
# LO SOLUTION( 101)      727.9505731659
# LO SOLUTION( 501)      6846.330143698
# LO SOLUTION(1001)      17176.03828316
# LO SOLUTION(5001)      
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
        pb.pbclass = "C-OOR2-AN-V-V"
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

    elseif action == "eEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EZ = exp(EV_[1])
        f_   = EZ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EZ
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = EZ
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

