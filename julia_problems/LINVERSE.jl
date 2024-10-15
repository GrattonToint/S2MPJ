function LINVERSE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
#    Problem : LINVERSE
#    *********
# 
#    The problem is to find the positive definite lower bidiagonal
#    matrix L such that the matrix L(inv)L(inv-transp) best approximates,
#    in the Frobenius norm, a given symmetric target matrix T.
#    More precisely, one is  interested in the positive definite lower
#    bidiagonal L such that
# 
#         || L T L(transp) - I ||     is minimum.
#                                F
# 
#    The positive definite character of L is imposed by requiring
#    that all its diagonal entries to be at least equal to EPSILON,
#    a strictly positive real number.
# 
#    Many variants of the problem can be obtained by varying the target
#    matrix T and the scalar EPSILON.  In the present problem,
#    a) T is chosen to be pentadiagonal with T(i,j) = sin(i)cos(j) (j .leq. i)
#    b) EPSILON = 1.D-8
# 
#    Source:
#    Ph. Toint, private communication, 1991.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "C-SBR2-AN-V-0"
# 
#    Dimension of the matrix
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER  n = 19    original value
# IE N                   100            $-PARAMETER  n = 199
# IE N                   500            $-PARAMETER  n = 999
# IE N                   1000           $-PARAMETER  n = 1999
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LINVERSE"

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
        v_["EPSILON"] = 1.0e-8
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["N-1"] = -1+v_["N"]
        v_["N-2"] = -2+v_["N"]
        for J = Int64(v_["1"]):Int64(v_["N-2"])
            v_["J+2"] = 2+J
            v_["RJ"] = Float64(J)
            v_["COSJ"] = cos(v_["RJ"])
            for I = Int64(J):Int64(v_["J+2"])
                v_["RI"] = Float64(I)
                v_["SINI"] = sin(v_["RI"])
                v_["T"*string(I)*","*string(J)] = v_["SINI"]*v_["COSJ"]
            end
        end
        v_["RN-1"] = Float64(v_["N-1"])
        v_["SINI"] = sin(v_["RN-1"])
        v_["COSJ"] = cos(v_["RN-1"])
        v_["T"*string(Int64(v_["N-1"]))*","*string(Int64(v_["N-1"]))]  = (
              v_["SINI"]*v_["COSJ"])
        v_["RN"] = Float64(v_["N"])
        v_["SINI"] = sin(v_["RN"])
        v_["T"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"]))]  = (
              v_["SINI"]*v_["COSJ"])
        v_["COSJ"] = cos(v_["RN"])
        v_["T"*string(Int64(v_["N"]))*","*string(Int64(v_["N"]))]  = (
              v_["SINI"]*v_["COSJ"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("A"*string(I),ix_)
            arrset(pb.xnames,iv,"A"*string(I))
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("A"*string(Int64(v_["N"])),ix_)
        arrset(pb.xnames,iv,"A"*string(Int64(v_["N"])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for J = Int64(v_["1"]):Int64(v_["N-2"])
            v_["J+1"] = 1+J
            v_["J+2"] = 2+J
            ig,ig_,_ = s2mpj_ii("O"*string(J)*","*string(J),ig_)
            arrset(gtype,ig,"<>")
            ig,ig_,_ = s2mpj_ii("O"*string(Int64(v_["J+1"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(0.5))
            ig,ig_,_ = s2mpj_ii("O"*string(Int64(v_["J+2"]))*","*string(J),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(0.5))
        end
        ig,ig_,_  = (
              s2mpj_ii("O"*string(Int64(v_["N-1"]))*","*string(Int64(v_["N-1"])),ig_))
        arrset(gtype,ig,"<>")
        ig,ig_,_  = (
              s2mpj_ii("O"*string(Int64(v_["N"]))*","*string(Int64(v_["N-1"])),ig_))
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(0.5))
        ig,ig_,_  = (
              s2mpj_ii("O"*string(Int64(v_["N"]))*","*string(Int64(v_["N"])),ig_))
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["O"*string(I)*","*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["A"*string(I)]] = v_["EPSILON"]
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(-1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "S"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "S"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "U"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "A"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "V"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "W"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "B"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["4"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["I-2"] = -2+I
            ename = "S"*string(I)*","*string(Int64(v_["I-2"]))
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
            end
            vname = "A"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "S"*string(I)*","*string(Int64(v_["I-2"]))
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
            end
            vname = "A"*string(Int64(v_["I-2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "V"*string(I)*","*string(Int64(v_["I-2"]))
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
            end
            vname = "B"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "V"*string(I)*","*string(Int64(v_["I-2"]))
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
            end
            vname = "A"*string(Int64(v_["I-2"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            for J = Int64(v_["I-1"]):Int64(I)
                v_["J-1"] = -1+J
                ename = "S"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                end
                vname = "A"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "A"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "U"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                end
                vname = "A"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "B"*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "V"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                end
                vname = "B"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "A"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "W"*string(I)*","*string(J)
                ie,ie_,newelt = s2mpj_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"en2PR")
                    arrset(ielftype,ie,iet_["en2PR"])
                end
                vname = "B"*string(Int64(v_["I-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "B"*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(-1.0))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
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
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        ig = ig_["O"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
        ig = ig_["O"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["V"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
        ig = ig_["O"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["V"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
        ig = ig_["O"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["U"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["V"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["W"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
        ig = ig_["O"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["U"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["1"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["V"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["W"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["1"]))]))
        ig = ig_["O"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["U"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["V"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["3"]))*","*string(Int64(v_["2"]))]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["W"*string(Int64(v_["3"]))*","*string(Int64(v_["3"]))])
        loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["2"]))*","*string(Int64(v_["2"]))]))
        for I = Int64(v_["4"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["I-2"] = -2+I
            ig = ig_["O"*string(I)*","*string(Int64(v_["I-2"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(Int64(v_["I-2"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(Int64(v_["I-2"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["V"*string(I)*","*string(Int64(v_["I-2"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-2"]))]))
            ig = ig_["O"*string(I)*","*string(Int64(v_["I-1"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(Int64(v_["I-1"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(Int64(v_["I-1"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U"*string(I)*","*string(Int64(v_["I-1"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(Int64(v_["I-2"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["V"*string(I)*","*string(Int64(v_["I-1"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-1"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["W"*string(I)*","*string(Int64(v_["I-1"]))])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-2"]))]))
            ig = ig_["O"*string(I)*","*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)*","*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(I)]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U"*string(I)*","*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(Int64(v_["I-1"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["V"*string(I)*","*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(I)*","*string(Int64(v_["I-1"]))]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["W"*string(I)*","*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["T"*string(Int64(v_["I-1"]))*","*string(Int64(v_["I-1"]))]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           6.00000000
# LO SOLTN(100)          68.0000000
# LO SOLTN(500)          340.000000
# LO SOLTN(1000)         ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-AN-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

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

