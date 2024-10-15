function JUNKTURN(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : JUNKTURN
#    *********
# 
#    The spacecraft orientation problem by Junkins and Turner. This is a
#    nonlinear optimal control problem.
# 
#    The problem is not convex.
# 
#    Source:
#    A.I Tyatushkin, A.I. Zholudev and N. M. Erinchek,
#    "The gradient method for solving optimal control problems with phase
#    constraints", 
#    in "System Modelling and Optimization", P. Kall, ed., pp. 456--464,
#    Springer Verlag, Lecture Notes in Control and Information Sciences 180, 1992.
#    This reference itself refers to:
#    I.L. Junkins and I.D. Turner,
#    "Optimal continuous torque attitude maneuvers",
#    AIAA/AAS Astrodynamics Conference, Palo Alto, 1978.
# 
#    SIF input: Ph. Toint, February 1994.
# 
#    classification = "C-QQR2-MN-V-V"
# 
#    Number of discretized points in [0,100] - 1.
#    The number of variables is    10 * ( N + 1 )
#    The number of constraints is  7 * N
#    N should be large enough to ensure feasibility.
# 
#       Alternative values for the SIF file parameters:
# IE N                   50             $-PARAMETER n =     510, m =    350
# IE N                   100            $-PARAMETER n =    1010, m =    700
# IE N                   500            $-PARAMETER n =    5010, m =   3500
# IE N                   1000           $-PARAMETER n =   10010, m =   7000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "JUNKTURN"

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
            v_["N"] = Int64(5);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   10000          $-PARAMETER n =  100010, m =  70000
# IE N                   20000          $-PARAMETER n =  200010, m = 140000
# IE N                   100000         $-PARAMETER n = 1000010, m = 700000
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = 100.0/v_["RN"]
        v_["H/2"] = 0.5*v_["H"]
        v_["6H/5"] = 1.2*v_["H"]
        v_["SH"] = 1.0909*v_["H"]
        v_["S1H"] = -0.08333*v_["H"]
        v_["S2H"] = 0.18182*v_["H"]
        v_["-H/2"] = -0.5*v_["H"]
        v_["-H"] = -1.0*v_["H"]
        v_["-H/10"] = -0.1*v_["H"]
        v_["H/4"] = 0.25*v_["H"]
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
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(T),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(T))
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
        for T = Int64(v_["1"]):Int64(v_["N"])
            v_["T-1"] = -1+T
            for I = Int64(v_["1"]):Int64(v_["7"])
                ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(T),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"C"*string(I)*","*string(T))
                iv = ix_["X"*string(I)*","*string(T)]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(I)*","*string(Int64(v_["T-1"]))]
                pbm.A[ig,iv] += Float64(1.0)
            end
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["5"]))*","*string(T),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["5"]))*","*string(T))
            iv = ix_["U"*string(Int64(v_["1"]))*","*string(T)]
            pbm.A[ig,iv] += Float64(v_["H"])
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["6"]))*","*string(T),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["6"]))*","*string(T))
            iv = ix_["U"*string(Int64(v_["2"]))*","*string(T)]
            pbm.A[ig,iv] += Float64(v_["6H/5"])
            ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["7"]))*","*string(T),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(Int64(v_["7"]))*","*string(T))
            iv = ix_["U"*string(Int64(v_["3"]))*","*string(T)]
            pbm.A[ig,iv] += Float64(v_["SH"])
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
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]] = 1.0
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]] = 1.0
        pb.xlower[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["3"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["3"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["4"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["4"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["5"]))*","*string(Int64(v_["0"]))]] = 0.01
        pb.xupper[ix_["X"*string(Int64(v_["5"]))*","*string(Int64(v_["0"]))]] = 0.01
        pb.xlower[ix_["X"*string(Int64(v_["6"]))*","*string(Int64(v_["0"]))]]  = (
              0.005)
        pb.xupper[ix_["X"*string(Int64(v_["6"]))*","*string(Int64(v_["0"]))]]  = (
              0.005)
        pb.xlower[ix_["X"*string(Int64(v_["7"]))*","*string(Int64(v_["0"]))]]  = (
              0.001)
        pb.xupper[ix_["X"*string(Int64(v_["7"]))*","*string(Int64(v_["0"]))]]  = (
              0.001)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["N"]))]]  = (
              0.43047)
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["N"]))]]  = (
              0.43047)
        pb.xlower[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["N"]))]]  = (
              0.70106)
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["N"]))]]  = (
              0.70106)
        pb.xlower[ix_["X"*string(Int64(v_["3"]))*","*string(Int64(v_["N"]))]]  = (
              0.0923)
        pb.xupper[ix_["X"*string(Int64(v_["3"]))*","*string(Int64(v_["N"]))]]  = (
              0.0923)
        pb.xlower[ix_["X"*string(Int64(v_["4"]))*","*string(Int64(v_["N"]))]]  = (
              0.56098)
        pb.xupper[ix_["X"*string(Int64(v_["4"]))*","*string(Int64(v_["N"]))]]  = (
              0.56098)
        pb.xlower[ix_["X"*string(Int64(v_["5"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["5"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["6"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["6"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["7"]))*","*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["7"]))*","*string(Int64(v_["N"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
              Float64(1.0))
        pb.x0[ix_["X"*string(Int64(v_["5"]))*","*string(Int64(v_["0"]))]]  = (
              Float64(0.01))
        pb.x0[ix_["X"*string(Int64(v_["6"]))*","*string(Int64(v_["0"]))]]  = (
              Float64(0.005))
        pb.x0[ix_["X"*string(Int64(v_["7"]))*","*string(Int64(v_["0"]))]]  = (
              Float64(0.001))
        pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["N"]))]]  = (
              Float64(0.43047))
        pb.x0[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["N"]))]]  = (
              Float64(0.70106))
        pb.x0[ix_["X"*string(Int64(v_["3"]))*","*string(Int64(v_["N"]))]]  = (
              Float64(0.0923))
        pb.x0[ix_["X"*string(Int64(v_["4"]))*","*string(Int64(v_["N"]))]]  = (
              Float64(0.56098))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for T = Int64(v_["0"]):Int64(v_["N"])
            ename = "U1S"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(Int64(v_["1"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "U2S"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(Int64(v_["2"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "U3S"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(Int64(v_["3"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for T = Int64(v_["1"]):Int64(v_["N"])
            ename = "P15"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["1"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P16"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["1"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P17"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["1"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P25"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["2"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P26"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["2"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P27"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["2"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P35"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["3"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P36"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["3"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P37"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["3"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P45"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["4"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P46"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["4"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P47"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["4"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P56"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P57"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["5"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "P67"*string(T)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(Int64(v_["6"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["7"]))*","*string(T)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["U1S"*string(Int64(v_["0"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["H/4"]))
        for T = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U1S"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U2S"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U3S"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
        end
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["U1S"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["H/4"]))
        for T = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["C"*string(Int64(v_["1"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P25"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P36"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P47"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            ig = ig_["C"*string(Int64(v_["2"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P15"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P37"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P46"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            ig = ig_["C"*string(Int64(v_["3"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P16"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P27"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P45"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            ig = ig_["C"*string(Int64(v_["4"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P17"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P26"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P35"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            ig = ig_["C"*string(Int64(v_["5"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P67"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["S1H"]))
            ig = ig_["C"*string(Int64(v_["6"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P57"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/10"]))
            ig = ig_["C"*string(Int64(v_["7"]))*","*string(T)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P56"*string(T)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["S2H"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(500)          7.417771100D-5
# LO SOLTN(1000)         1.224842784D-5
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
        pb.pbclass = "C-QQR2-MN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

