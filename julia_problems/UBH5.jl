function UBH5(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : UBH5
#    *********
# 
#    The problem is to minimize the integral of the control magnitude needed
#    to bring a vehicle, from given position and velocity, to the origin with
#    zero velocity in a fixed amount of time.  The controls are the components
#    of the vehicle acceleration. The discretization uses the trapezoidal rule.
#    This version of the problem is a variant of UBH1, where the cumulative
#    value of the objective is maintained as an additional state variable.
# 
#    The problem is convex.
# 
#    Source: unscaled problem 5 
#    (ODE = 1, CLS = 2, GRD = 1, MET = T, SEED = 0.) in
#    J.T. Betts and W.P. Huffman,
#    "Sparse Nonlinear Programming Test Problems (Release 1.0)",
#    Boeing Computer services, Seattle, July 1993.
# 
#    SIF input: Ph.L. Toint, October 1993.
# 
#    classification = "C-LQR2-MN-V-V"
# 
#    Number of grid points
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n=100, m=70    original value
# IE N                   100            $-PARAMETER n=1000, m=700
# IE N                   500            $-PARAMETER n=5000, m=3500
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "UBH5"

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
# IE N                   1000           $-PARAMETER n=10000, m=7000
# IE N                   2000           $-PARAMETER n=20000, m=14000
        v_["T0"] = 0.0
        v_["TF"] = 1000.0
        v_["RN"] = Float64(v_["N"])
        v_["TTIME"] = v_["TF"]-v_["T0"]
        v_["K"] = v_["TTIME"]/v_["RN"]
        v_["-K/2"] = -0.5*v_["K"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["7"] = 7
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["7"])
            for T = Int64(v_["0"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("Y"*string(I)*","*string(T),ix_)
                arrset(pb.xnames,iv,"Y"*string(I)*","*string(T))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["3"])
            for T = Int64(v_["0"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("U"*string(I)*","*string(T),ix_)
                arrset(pb.xnames,iv,"U"*string(I)*","*string(T))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Y"*string(Int64(v_["7"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["3"])
            v_["I+3"] = 3+I
            for T = Int64(v_["1"]):Int64(v_["N"])
                v_["T-1"] = -1+T
                ig,ig_,_ = s2mpj_ii("S"*string(I)*","*string(T),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"S"*string(I)*","*string(T))
                iv = ix_["Y"*string(I)*","*string(T)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(I)*","*string(Int64(v_["T-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["Y"*string(Int64(v_["I+3"]))*","*string(Int64(v_["T-1"]))]
                pbm.A[ig,iv] += Float64(v_["-K/2"])
                iv = ix_["Y"*string(Int64(v_["I+3"]))*","*string(T)]
                pbm.A[ig,iv] += Float64(v_["-K/2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["3"])
            v_["I+3"] = 3+I
            for T = Int64(v_["1"]):Int64(v_["N"])
                v_["T-1"] = -1+T
                ig,ig_,_ = s2mpj_ii("S"*string(Int64(v_["I+3"]))*","*string(T),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"S"*string(Int64(v_["I+3"]))*","*string(T))
                iv = ix_["Y"*string(Int64(v_["I+3"]))*","*string(T)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(Int64(v_["I+3"]))*","*string(Int64(v_["T-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                ig,ig_,_ = s2mpj_ii("S"*string(Int64(v_["I+3"]))*","*string(T),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"S"*string(Int64(v_["I+3"]))*","*string(T))
                iv = ix_["U"*string(I)*","*string(Int64(v_["T-1"]))]
                pbm.A[ig,iv] += Float64(v_["-K/2"])
                ig,ig_,_ = s2mpj_ii("S"*string(Int64(v_["I+3"]))*","*string(T),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"S"*string(Int64(v_["I+3"]))*","*string(T))
                iv = ix_["U"*string(I)*","*string(T)]
                pbm.A[ig,iv] += Float64(v_["-K/2"])
            end
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            v_["T-1"] = -1+T
            ig,ig_,_ = s2mpj_ii("S"*string(Int64(v_["7"]))*","*string(T),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"S"*string(Int64(v_["7"]))*","*string(T))
            iv = ix_["Y"*string(Int64(v_["7"]))*","*string(T)]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["Y"*string(Int64(v_["7"]))*","*string(Int64(v_["T-1"]))]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["3"])
            for T = Int64(v_["0"]):Int64(v_["N"])
                pb.xlower[ix_["U"*string(I)*","*string(T)]] = -1.0
                pb.xupper[ix_["U"*string(I)*","*string(T)]] = 1.0
            end
        end
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xlower[ix_["Y"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xupper[ix_["Y"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xlower[ix_["Y"*string(Int64(v_["3"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xupper[ix_["Y"*string(Int64(v_["3"]))*","*string(Int64(v_["0"]))]]  = (
              1000.0)
        pb.xlower[ix_["Y"*string(Int64(v_["4"]))*","*string(Int64(v_["0"]))]] = (-
             10.0)
        pb.xupper[ix_["Y"*string(Int64(v_["4"]))*","*string(Int64(v_["0"]))]] = (-
             10.0)
        pb.xlower[ix_["Y"*string(Int64(v_["5"]))*","*string(Int64(v_["0"]))]] = 10.0
        pb.xupper[ix_["Y"*string(Int64(v_["5"]))*","*string(Int64(v_["0"]))]] = 10.0
        pb.xlower[ix_["Y"*string(Int64(v_["6"]))*","*string(Int64(v_["0"]))]] = (-
             10.0)
        pb.xupper[ix_["Y"*string(Int64(v_["6"]))*","*string(Int64(v_["0"]))]] = (-
             10.0)
        pb.xlower[ix_["Y"*string(Int64(v_["7"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["7"]))*","*string(Int64(v_["0"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["6"])
            pb.xlower[ix_["Y"*string(I)*","*string(Int64(v_["N"]))]] = 0.0
            pb.xupper[ix_["Y"*string(I)*","*string(Int64(v_["N"]))]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["0"]):Int64(v_["N"])
            for I = Int64(v_["1"]):Int64(v_["3"])
                ename = "E"*string(I)*","*string(T)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "U"*string(I)*","*string(T)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for T = Int64(v_["1"]):Int64(v_["N"])
            v_["T-1"] = -1+T
            for I = Int64(v_["1"]):Int64(v_["3"])
                ig = ig_["S"*string(Int64(v_["7"]))*","*string(T)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(Int64(v_["T-1"]))])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-K/2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(T)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-K/2"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(10)           1.14735202967
# LO SOLTN(100)          1.11631518169
# LO SOLTN(1000)         1.11598643493
# LO SOLTN(2000)         1.11587382445
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
        pb.pbclass = "C-LQR2-MN-V-V"
        pb.x0          = zeros(Float64,pb.n)
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

