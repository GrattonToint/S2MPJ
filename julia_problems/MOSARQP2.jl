function MOSARQP2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MOSARQP2
#    *********
# 
#    A convex quadratic problem, with variable dimensions.
#    In this problem, a third of the linear constraints are active at the
#    solution. 
# 
#    Source:
#    J.L. Morales-Perez and R.W.H. Sargent,
#    "On the implementation and performance of an interior point method for
#    large sparse convex quadratic programming",
#    Centre for Process Systems Engineering, Imperial College, London,
#    November 1991.
# 
#    SIF input: Ph. Toint, August 1993.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-QLR2-AN-V-V"
# 
#    Problem variants: these are distinguished by the triplet ( N, M, COND ),
#    where: - N (nb of variables) must be a multiple of 3
#             and have an integer square root
#           - M (nb of constraints) must be at least sqrt(N) 
#             and at most N - sqrt(N)
#           - COND (problem conditioning) is a positive integer
#    Except for the first, the instances suggested are those used by Morales
#    and Sargent.
# 
#       Alternative values for the SIF file parameters:
# IE N                   36             $-PARAMETER
# IE M                   10             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   100            $-PARAMETER     original value
# IE M                   10             $-PARAMETER     original value
# RE COND                3.0            $-PARAMETER     original value
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   30             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   60             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   90             $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   120            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   300            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   900            $-PARAMETER
# IE M                   600            $-PARAMETER
# RE COND                3.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                1.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# IE M                   700            $-PARAMETER
# RE COND                2.0            $-PARAMETER
# 
# IE N                   2500           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MOSARQP2"

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
            v_["N"] = Int64(36);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE M                   700            $-PARAMETER
        if nargin<2
            v_["M"] = Int64(10);  #  SIF file default value
        else
            v_["M"] = Int64(args[2]);
        end
# RE COND                3.0            $-PARAMETER
        if nargin<3
            v_["COND"] = Float64(2.0);  #  SIF file default value
        else
            v_["COND"] = Float64(args[3]);
        end
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["N-1"] = -1+v_["N"]
        v_["N-2"] = -2+v_["N"]
        v_["RN-1"] = Float64(v_["N-1"])
        v_["RN"] = Float64(v_["N"])
        v_["M-1"] = -1+v_["M"]
        v_["RNP"] = 0.1+v_["RN"]
        v_["RRTN"] = sqrt(v_["RNP"])
        v_["RTN"] = trunc(Int,v_["RRTN"])
        v_["RTN-1"] = -1+v_["RTN"]
        v_["RTN-2"] = -2+v_["RTN"]
        v_["RTN+1"] = 1+v_["RTN"]
        v_["2RTN-1"] = v_["RTN"]+v_["RTN-1"]
        v_["2RTN"] = v_["RTN"]+v_["RTN"]
        v_["M-RTN+1"] = v_["M"]-v_["RTN-1"]
        for I = Int64(v_["1"]):Int64(v_["3"]):Int64(v_["N-2"])
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            v_["XC"*string(I)] = -1.0
            v_["XC"*string(Int64(v_["I+1"]))] = -1.0
            v_["XC"*string(Int64(v_["I+2"]))] = 1.0
        end
        v_["XC"*string(Int64(v_["N-1"]))] = -1.0
        v_["XC"*string(Int64(v_["N"]))] = 1.0
        v_["NNZ"] = 10
        v_["Y1"] = -0.3569732
        v_["Y2"] = 0.9871576
        v_["Y3"] = 0.5619363
        v_["Y4"] = -0.1984624
        v_["Y5"] = 0.4653328
        v_["Y6"] = 0.7364367
        v_["Y7"] = -0.4560378
        v_["Y8"] = -0.6457813
        v_["Y9"] = -0.0601357
        v_["Y10"] = 0.1035624
        v_["NZ1"] = 0.68971452
        v_["NZ2"] = 0.13452678
        v_["NZ3"] = 0.51234678
        v_["NZ4"] = 0.76591423
        v_["NZ5"] = 0.20857854
        v_["NZ6"] = 0.85672348
        v_["NZ7"] = 0.04356789
        v_["NZ8"] = 0.44692743
        v_["NZ9"] = 0.30136413
        v_["NZ10"] = 0.91367489
        v_["YN2"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["NNZ"])
            v_["RKI"] = v_["NZ"*string(I)]*v_["RN"]
            v_["K"*string(I)] = 1.1+v_["RKI"]
            v_["TMP"] = v_["Y"*string(I)]*v_["Y"*string(I)]
            v_["YN2"] = v_["YN2"]+v_["TMP"]
        end
        v_["-2/YN2"] = -2.0/v_["YN2"]
        v_["4/YN4"] = v_["-2/YN2"]*v_["-2/YN2"]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TMP"] = v_["RI-1"]/v_["RN-1"]
            v_["TMP"] = v_["TMP"]*v_["COND"]
            v_["D"*string(I)] = exp(v_["TMP"])
        end
        v_["YDY"] = 0.0
        v_["YXC"] = 0.0
        v_["YDXC"] = 0.0
        for I = Int64(v_["1"]):Int64(v_["NNZ"])
            v_["RKI"] = v_["K"*string(I)]
            v_["KI"] = trunc(Int,v_["RKI"])
            v_["DY"*string(I)] = v_["Y"*string(I)]*v_["D"*string(Int64(v_["KI"]))]
            v_["TMP"] = v_["DY"*string(I)]*v_["Y"*string(I)]
            v_["YDY"] = v_["YDY"]+v_["TMP"]
            v_["TMP"] = v_["Y"*string(I)]*v_["XC"*string(Int64(v_["KI"]))]
            v_["YXC"] = v_["YXC"]+v_["TMP"]
            v_["TMP"] = v_["DY"*string(I)]*v_["XC"*string(Int64(v_["KI"]))]
            v_["YDXC"] = v_["YDXC"]+v_["TMP"]
        end
        v_["AA"] = v_["-2/YN2"]*v_["YXC"]
        v_["DD"] = v_["4/YN4"]*v_["YDY"]
        v_["BB"] = v_["DD"]*v_["YXC"]
        v_["CC"] = v_["-2/YN2"]*v_["YDXC"]
        v_["BB+CC"] = v_["BB"]+v_["CC"]
        v_["DD/2"] = 0.5*v_["DD"]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["C"*string(I)] = v_["D"*string(I)]*v_["XC"*string(I)]
        end
        for I = Int64(v_["1"]):Int64(v_["NNZ"])
            v_["RKI"] = v_["K"*string(I)]
            v_["KI"] = trunc(Int,v_["RKI"])
            v_["TMP"] = v_["DY"*string(I)]*v_["AA"]
            v_["C"*string(Int64(v_["KI"]))] = v_["C"*string(Int64(v_["KI"]))]+v_["TMP"]
            v_["TMP"] = v_["Y"*string(I)]*v_["BB+CC"]
            v_["C"*string(Int64(v_["KI"]))] = v_["C"*string(Int64(v_["KI"]))]+v_["TMP"]
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
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["C"*string(I)])
        end
        ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CS"*string(Int64(v_["1"])))
        iv = ix_["X"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(4.0)
        ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CS"*string(Int64(v_["1"])))
        iv = ix_["X"*string(Int64(v_["RTN+1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        for I = Int64(v_["2"]):Int64(v_["RTN-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["I+RTN"] = I+v_["RTN"]
            ig,ig_,_ = s2mpj_ii("CS"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(4.0)
            iv = ix_["X"*string(Int64(v_["I+RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["RTN"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CS"*string(Int64(v_["RTN"])))
        iv = ix_["X"*string(Int64(v_["RTN"]))]
        pbm.A[ig,iv] += Float64(4.0)
        ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["RTN"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CS"*string(Int64(v_["RTN"])))
        iv = ix_["X"*string(Int64(v_["RTN-1"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X"*string(Int64(v_["2RTN"]))]
        pbm.A[ig,iv] += Float64(-1.0)
        v_["JS"] = v_["RTN"]
        for J = Int64(v_["RTN+1"]):Int64(v_["RTN"]):Int64(v_["M-RTN+1"])
            v_["J+1"] = 1+J
            v_["JS"] = J+v_["RTN-1"]
            v_["JS-1"] = -1+v_["JS"]
            v_["J-RTN"] = J-v_["RTN"]
            v_["J+RTN"] = J+v_["RTN"]
            ig,ig_,_ = s2mpj_ii("CS"*string(J),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(J))
            iv = ix_["X"*string(J)]
            pbm.A[ig,iv] += Float64(4.0)
            iv = ix_["X"*string(Int64(v_["J+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["J-RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["J+RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            for I = Int64(v_["J+1"]):Int64(v_["JS-1"])
                v_["I+1"] = 1+I
                v_["I-1"] = -1+I
                v_["I+RTN"] = I+v_["RTN"]
                v_["I-RTN"] = I-v_["RTN"]
                ig,ig_,_ = s2mpj_ii("CS"*string(I),ig_)
                arrset(gtype,ig,">=")
                arrset(pb.cnames,ig,"CS"*string(I))
                iv = ix_["X"*string(I)]
                pbm.A[ig,iv] += Float64(4.0)
                iv = ix_["X"*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(Int64(v_["I+1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(Int64(v_["I-RTN"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["X"*string(Int64(v_["I+RTN"]))]
                pbm.A[ig,iv] += Float64(-1.0)
            end
            v_["JS+RTN"] = v_["JS"]+v_["RTN"]
            v_["JS-RTN"] = v_["JS"]-v_["RTN"]
            ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["JS"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(Int64(v_["JS"])))
            iv = ix_["X"*string(Int64(v_["JS"]))]
            pbm.A[ig,iv] += Float64(4.0)
            iv = ix_["X"*string(Int64(v_["JS-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["JS"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(Int64(v_["JS"])))
            iv = ix_["X"*string(Int64(v_["JS-RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["JS+RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        v_["K"] = 1+v_["JS"]
        for I = Int64(v_["K"]):Int64(v_["M"]):Int64(v_["M"])
            v_["K+1"] = 1+v_["K"]
            v_["K+RTN"] = v_["K"]+v_["RTN"]
            v_["K-RTN"] = v_["K"]-v_["RTN"]
            ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["K"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(Int64(v_["K"])))
            iv = ix_["X"*string(Int64(v_["K"]))]
            pbm.A[ig,iv] += Float64(4.0)
            iv = ix_["X"*string(Int64(v_["K+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("CS"*string(Int64(v_["K"])),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(Int64(v_["K"])))
            iv = ix_["X"*string(Int64(v_["K-RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["K+RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        v_["K"] = 1+v_["K"]
        for I = Int64(v_["K"]):Int64(v_["M"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["I+RTN"] = I+v_["RTN"]
            v_["I-RTN"] = I-v_["RTN"]
            ig,ig_,_ = s2mpj_ii("CS"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CS"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(4.0)
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["I-RTN"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["X"*string(Int64(v_["I+RTN"]))]
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
        pbm.gconst[ig_["CS"*string(Int64(v_["1"]))]] = Float64(0.5)
        pbm.gconst[ig_["CS"*string(Int64(v_["RTN"]))]] = Float64(0.5)
        v_["K"] = v_["RTN+1"]
        for J = Int64(v_["RTN+1"]):Int64(v_["RTN"]):Int64(v_["M-RTN+1"])
            v_["K"] = 1+v_["K"]
            for I = Int64(v_["1"]):Int64(v_["RTN-2"])
                v_["K"] = 1+v_["K"]
                pbm.gconst[ig_["CS"*string(Int64(v_["K"]))]] = Float64(-0.5)
            end
            v_["K"] = 1+v_["K"]
        end
        v_["K"] = 1+v_["K"]
        for J = Int64(v_["K"]):Int64(v_["M"])
            pbm.gconst[ig_["CS"*string(J)]] = Float64(-0.5)
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "XSQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["NNZ"])
            v_["RKI"] = v_["K"*string(I)]
            v_["KI"] = trunc(Int,v_["RKI"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["RKJ"] = v_["K"*string(J)]
                v_["KJ"] = trunc(Int,v_["RKJ"])
                ename = "P"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
                vname = "X"*string(Int64(v_["KI"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["KJ"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["TMP"] = 0.5*v_["D"*string(I)]
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["TMP"]))
        end
        for I = Int64(v_["1"]):Int64(v_["NNZ"])
            v_["RKI"] = v_["K"*string(I)]
            v_["KI"] = trunc(Int,v_["RKI"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["TMP"] = v_["DY"*string(I)]*v_["Y"*string(J)]
                v_["WIJ"] = v_["TMP"]*v_["-2/YN2"]
                v_["TMP"] = v_["DY"*string(J)]*v_["Y"*string(I)]
                v_["TMP"] = v_["TMP"]*v_["-2/YN2"]
                v_["WIJ"] = v_["WIJ"]+v_["TMP"]
                v_["TMP"] = v_["Y"*string(I)]*v_["Y"*string(J)]
                v_["TMP"] = v_["TMP"]*v_["DD"]
                v_["WIJ"] = v_["WIJ"]+v_["TMP"]
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["P"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["WIJ"]))
            end
            v_["TMP"] = v_["DY"*string(I)]*v_["Y"*string(I)]
            v_["WII"] = v_["TMP"]*v_["-2/YN2"]
            v_["TMP"] = v_["Y"*string(I)]*v_["Y"*string(I)]
            v_["TMP"] = v_["TMP"]*v_["DD/2"]
            v_["WII"] = v_["WII"]+v_["TMP"]
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["KI"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["WII"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(  36, 10,2)   -35.69811798
# LO SOLTN( 900, 30,1)   -509.8245900
# LO SOLTN( 900, 30,2)   -950.8404853
# LO SOLTN( 900, 30,3)   -1896.596722
# LO SOLTN( 900, 60,1)   -504.3600140
# LO SOLTN( 900, 60,2)   -945.1134463
# LO SOLTN( 900, 60,3)   -1890.602184
# LO SOLTN( 900, 90,1)   -498.9518964
# LO SOLTN( 900, 90,2)   -939.2704526
# LO SOLTN( 900, 90,3)   -1884.291256
# LO SOLTN( 900,120,1)   -493.5058050
# LO SOLTN( 900,120,2)   -933.1963138
# LO SOLTN( 900,120,3)   -1877.513644
# LO SOLTN( 900,300,1)   -457.1185630
# LO SOLTN( 900,300,2)   -887.3869230
# LO SOLTN( 900,300,3)   -1819.655008
# LO SOLTN( 900,600,1)   -377.5813314
# LO SOLTN( 900,600,2)   -755.0919955
# LO SOLTN( 900,600,3)   -1597.482277
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-QLR2-AN-V-V"
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

