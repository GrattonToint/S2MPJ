function NUFFIELD(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : NUFFIELD
#    *********
# 
#    A problem from economics.
#    Maximize a 2-D integral representing consumer surplus subject to 
#    linear and quadratic constraints representing incentive compatibility
# 
#    Let v( . , . ) : R^2 -> R, Omega = [a,a+1] x [a,a+1], and
#    the corners A, B, C, D be as follows:
# 
#            (a+1,a+1)
#        A *-----* B
#          |     |
#          |     |
#        D *-----* C
#        (a,a)  
# 
#    The problem is to maximize
# 
#       (a+1) line integral_{AB U BC} v(w)dw 
#        - a line integral_{CD U DA} v(w)dw
#        - 3 volume integral_{Omega} v(w)dw
# 
#    subject to v being symmetric (i.e., v(x,y) = v(y,x))
#               v(a,a) = 0
#               nabla_w v(w) >= 0
#               < e, nabla_w v(w) > <= 1
#         and   nabla_ww v(w) positive definite
# 
#    this last constraint is guaranteed by ensuring that
# 
#               d^2 v/dx^2 >= 0
#               d^2 v/dy^2 >= 0
#               ( d^2 v/dx^2 )( d^2 v/dy^2 ) >= ( d^2 v/dxdy )^2
# 
#    Symmetry is ensured by only considering v(x,y) for x <= y
# 
#    Here v(x,y) is the consumer surplus. that is if the consumer values good 
#    1 at x pounds and good 2 at y pounds then they will have a utility 
#    equivalent to v(x,y) pounds after being faced with the optimal monopoly 
#    pricing strategy. (Apparently, from this we can infer what the optimal 
#    pricing strategy was... ).
# 
#    More background is available from
# 
#    "Optimal Selling Strategies: When to haggle, when to hold firm",
#      Riley and Zeckhauser. The Quarterly Journal of Economics, 1983, and
# 
#    "Multidimensional Incentive Compatibility and Mechanism Design", 
#      McAfee and McMillan. The Journal of Economic Theory, 1988.
# 
#    Source: John Thanassoulis <john.thanassoulis@nuffield.oxford.ac.uk>
# 
#    Standard finite-differences are used to ap[proximate derivatives, and 
#    1- and 2-D trapezoidal rules to approximate integrals
# 
#    SIF input: Nick Gould, February 2001
# 
#    classification = "C-LQR2-AN-V-V"
# 
#    The parameter a
# 
#       Alternative values for the SIF file parameters:
# RE A                   5.0            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "NUFFIELD"

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
            v_["A"] = Float64(5.0);  #  SIF file default value
        else
            v_["A"] = Float64(args[1]);
        end
# IE N                   10            $-PARAMETER
# IE N                   20            $-PARAMETER
# IE N                   30            $-PARAMETER
# IE N                   40            $-PARAMETER
        if nargin<2
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[2]);
        end
# IE N                   100           $-PARAMETER
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = 1.0/v_["RN"]
        v_["1/H"] = v_["RN"]
        v_["-1/H"] = -1.0*v_["1/H"]
        v_["H**2"] = v_["H"]*v_["H"]
        v_["1/H**2"] = v_["1/H"]*v_["1/H"]
        v_["-2/H**2"] = -2.0*v_["1/H**2"]
        v_["1/H**4"] = v_["1/H**2"]*v_["1/H**2"]
        v_["A+1"] = 1.0+v_["A"]
        v_["-A-1"] = -1.0*v_["A+1"]
        v_["C2"] = 3.0*v_["H"]
        v_["C3"] = 0.5*v_["C2"]
        v_["C4"] = v_["C3"]+v_["A"]
        v_["C1"] = v_["C3"]+v_["-A-1"]
        v_["C5"] = -1.0+v_["C3"]
        v_["C5"] = 0.5*v_["C5"]
        v_["C6"] = 0.5*v_["C3"]
        v_["C6"] = v_["C6"]+v_["-A-1"]
        v_["C6"] = 0.5*v_["C6"]
        v_["C1"] = v_["C1"]*v_["H"]
        v_["C2"] = v_["C2"]*v_["H"]
        v_["C3"] = v_["C3"]*v_["H"]
        v_["C4"] = v_["C4"]*v_["H"]
        v_["C5"] = v_["C5"]*v_["H"]
        v_["C6"] = v_["C6"]*v_["H"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(I)
                iv,ix_,_ = s2mpj_ii("V"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"V"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["V"*string(Int64(v_["N"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["C1"])
        end
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("OBJ",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["C2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["V"*string(I)*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["C3"])
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["V"*string(I)*","*string(Int64(v_["0"]))]
            pbm.A[ig,iv] += Float64(v_["C4"])
        end
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(v_["C5"])
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["C6"])
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            for J = Int64(v_["0"]):Int64(I)
                ig,ig_,_ = s2mpj_ii("VX"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"VX"*string(I)*","*string(J))
                iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-1/H"])
                ig,ig_,_ = s2mpj_ii("VV"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"VV"*string(I)*","*string(J))
                iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-1/H"])
            end
        end
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("VX"*string(Int64(v_["N"]))*","*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"VX"*string(Int64(v_["N"]))*","*string(J))
            iv = ix_["V"*string(Int64(v_["N"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            ig,ig_,_ = s2mpj_ii("VX"*string(Int64(v_["N"]))*","*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"VX"*string(Int64(v_["N"]))*","*string(J))
            iv = ix_["V"*string(Int64(v_["N-1"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            ig,ig_,_ = s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(J),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(J))
            iv = ix_["V"*string(Int64(v_["N"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            ig,ig_,_ = s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(J),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(J))
            iv = ix_["V"*string(Int64(v_["N-1"]))*","*string(J)]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
        end
        ig,ig_,_  = (
              s2mpj_ii("VX"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"VX"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VX"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"VX"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(v_["-1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(v_["-1/H"])
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            for J = Int64(v_["0"]):Int64(v_["I-1"])
                v_["J+1"] = 1+J
                ig,ig_,_ = s2mpj_ii("VY"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"VY"*string(I)*","*string(J))
                iv = ix_["V"*string(I)*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(v_["1/H"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-1/H"])
                ig,ig_,_ = s2mpj_ii("VV"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"VV"*string(I)*","*string(J))
                iv = ix_["V"*string(I)*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(v_["1/H"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-1/H"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("VY"*string(I)*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"VY"*string(I)*","*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["V"*string(I)*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            ig,ig_,_ = s2mpj_ii("VV"*string(I)*","*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"VV"*string(I)*","*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["V"*string(I)*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
        end
        ig,ig_,_  = (
              s2mpj_ii("VY"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"VY"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VY"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"VY"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(v_["-1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["1/H"])
        ig,ig_,_  = (
              s2mpj_ii("VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"VV"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])))
        iv = ix_["V"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(v_["-1/H"])
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["0"]):Int64(v_["I-1"])
                ig,ig_,_ = s2mpj_ii("VXX"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"VXX"*string(I)*","*string(J))
                iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H**2"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2/H**2"])
                iv = ix_["V"*string(Int64(v_["I-1"]))*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["1/H**2"])
            end
        end
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["J-1"] = -1+J
                v_["J+1"] = 1+J
                ig,ig_,_ = s2mpj_ii("VYY"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"VYY"*string(I)*","*string(J))
                iv = ix_["V"*string(I)*","*string(Int64(v_["J+1"]))]
                pbm.A[ig,iv] += Float64(v_["1/H**2"])
                iv = ix_["V"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["-2/H**2"])
                iv = ix_["V"*string(I)*","*string(Int64(v_["J-1"]))]
                pbm.A[ig,iv] += Float64(v_["1/H**2"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("VXX"*string(I)*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"VXX"*string(I)*","*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H**2"])
            iv = ix_["V"*string(I)*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["-2/H**2"])
            iv = ix_["V"*string(I)*","*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["1/H**2"])
            ig,ig_,_ = s2mpj_ii("VYY"*string(I)*","*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"VYY"*string(I)*","*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H**2"])
            iv = ix_["V"*string(I)*","*string(I)]
            pbm.A[ig,iv] += Float64(v_["-2/H**2"])
            iv = ix_["V"*string(I)*","*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["1/H**2"])
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            for J = Int64(v_["1"]):Int64(I)
                ig,ig_,_ = s2mpj_ii("C"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"C"*string(I)*","*string(J))
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
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["0"]):Int64(I)
                pbm.gconst[ig_["VV"*string(I)*","*string(J)]] = Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["0"]))*","*string(Int64(v_["0"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eCONVEX", iet_)
        loaset(elftv,it,1,"VIP1J")
        loaset(elftv,it,2,"VIJP1")
        loaset(elftv,it,3,"VIJ")
        loaset(elftv,it,4,"VIM1J")
        loaset(elftv,it,5,"VIJM1")
        loaset(elftv,it,6,"VIPJP")
        loaset(elftv,it,7,"VIPJM")
        loaset(elftv,it,8,"VIMJM")
        loaset(elftv,it,9,"VIMJP")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(I)
                v_["J+1"] = 1+J
                v_["J-1"] = -1+J
                ename = "C"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eCONVEX")
                arrset(ielftype,ie,iet_["eCONVEX"])
                vname = "V"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIP1J",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIM1J",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIJP1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIJM1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIJ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIPJP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIMJM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIPJM",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "V"*string(Int64(v_["I-1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="VIMJP",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            for J = Int64(v_["1"]):Int64(I)
                ig = ig_["C"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["1/H**4"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solutions (may be local!)
# LO SOLTN               -2.512312500   $ (n=10)
# LO SOLTN               -2.512359371   $ (n=20)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-AN-V-V"
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

    elseif action == "eCONVEX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,9)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-2
        U_[1,4] = U_[1,4]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]-2
        U_[2,5] = U_[2,5]+1
        U_[3,6] = U_[3,6]+2.500000e-01
        U_[3,8] = U_[3,8]+2.500000e-01
        U_[3,9] = U_[3,9]-2.500000e-01
        U_[3,7] = U_[3,7]-2.500000e-01
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        f_   = IV_[1]*IV_[2]-IV_[3]*IV_[3]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[2]
            g_[2] = IV_[1]
            g_[3] = -2.0*IV_[3]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
                H_[3,3] = -2.0
                H_ = U_'*H_*U_
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

