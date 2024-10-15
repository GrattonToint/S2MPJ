function JNLBRNG1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : JNLBRNG1
#    *********
# 
#    The quadratic journal bearing problem (with excentricity = 0.1)
#    This is a variant of the problem stated in the report quoted below.
#    It corresponds to the problem as distributed in MINPACK-2.
# 
#    Source:
#    J. More' and G. Toraldo,
#    "On the Solution of Large Quadratic-Programming Problems with Bound
#    Constraints", 
#    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
# 
#    SIF input: Ph. Toint, Dec 1989.
#    modified by Peihuang Chen, according to MINPACK-2, Apr 1992
# 
#    classification = "C-QBR2-AY-V-0"
# 
#    The rectangle is discretized into (pt-1)(py-1) little rectangles. The
#    heights of the considered surface above the corners of these little
#    rectangles are the problem variables,  There are px*py of them.
# 
#    PT is the number of points along the T (\theta) side of the rectangle
#    PY is the number of points along the Y side of the rectangle
# 
#       Alternative values for the SIF file parameters:
# IE PT                  4              $-PARAMETER  n=16
# IE PY                  4              $-PARAMETER
# 
# IE PT                  10             $-PARAMETER  n=100
# IE PY                  10             $-PARAMETER
# 
# IE PT                  23             $-PARAMETER  n=529
# IE PY                  23             $-PARAMETER
# 
# IE PT                  32             $-PARAMETER  n=1024
# IE PY                  32             $-PARAMETER
# 
# IE PT                  34             $-PARAMETER  n=1156
# IE PY                  34             $-PARAMETER
# 
# IE PT                  75             $-PARAMETER  n=5625   original value
# IE PY                  75             $-PARAMETER           original value
# 
# IE PT                  100            $-PARAMETER  n=10000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "JNLBRNG1"

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
            v_["PT"] = Int64(5);  #  SIF file default value
        else
            v_["PT"] = Int64(args[1]);
        end
# IE PY                  100            $-PARAMETER
        if nargin<2
            v_["PY"] = Int64(5);  #  SIF file default value
        else
            v_["PY"] = Int64(args[2]);
        end
# IE PT                  125            $-PARAMETER  n=15625
# IE PY                  125            $-PARAMETER
# IE PT                  4              $-PARAMETER  n=16
# IE PY                  4              $-PARAMETER
        if nargin<3
            v_["EX"] = Float64(0.1);  #  SIF file default value
        else
            v_["EX"] = Float64(args[3]);
        end
        v_["PI/4"] = atan(1.0)
        v_["LT"] = 8.0*v_["PI/4"]
        v_["LY"] = 20.0
        v_["SIX"] = 6.0
        v_["PT-1"] = -1+v_["PT"]
        v_["RPT-1"] = Float64(v_["PT-1"])
        v_["HT1"] = 1.0/v_["RPT-1"]
        v_["HT"] = v_["HT1"]*v_["LT"]
        v_["1/HT"] = 1.0/v_["HT"]
        v_["PY-1"] = -1+v_["PY"]
        v_["RPY-1"] = Float64(v_["PY-1"])
        v_["HY1"] = 1.0/v_["RPY-1"]
        v_["HY"] = v_["HY1"]*v_["LY"]
        v_["1/HY"] = 1.0/v_["HY"]
        v_["HTHY"] = v_["HT"]*v_["HY"]
        v_["HT/HY"] = v_["HT"]*v_["1/HY"]
        v_["HY/HT"] = v_["HY"]*v_["1/HT"]
        v_["EXHTHY"] = v_["HTHY"]*v_["EX"]
        v_["CLINC"] = -1.0*v_["EXHTHY"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["PT"])
            for J = Int64(v_["1"]):Int64(v_["PY"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["PT-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HT"]
            v_["SXI1"] = sin(v_["XI1"])
            v_["COEFF"] = v_["SXI1"]*v_["CLINC"]
            for J = Int64(v_["2"]):Int64(v_["PY-1"])
                ig,ig_,_ = s2mpj_ii("G",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["COEFF"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["PT-1"])
            for J = Int64(v_["1"]):Int64(v_["PY-1"])
                ig,ig_,_ = s2mpj_ii("GR"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(2.0))
            end
        end
        for I = Int64(v_["2"]):Int64(v_["PT"])
            for J = Int64(v_["2"]):Int64(v_["PY"])
                ig,ig_,_ = s2mpj_ii("GL"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(2.0))
            end
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for J = Int64(v_["1"]):Int64(v_["PY"])
            pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xlower[ix_["X"*string(Int64(v_["PT"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["PT"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["PT-1"])
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["PY"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["PY"]))]] = 0.0
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["2"]):Int64(v_["PT-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HT"]
            v_["SXI1"] = sin(v_["XI1"])
            for J = Int64(v_["2"]):Int64(v_["PY-1"])
                pb.x0[ix_["X"*string(I)*","*string(J)]] = Float64(v_["SXI1"])
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["PT-1"])
            v_["I+1"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["PY-1"])
                v_["J+1"] = 1+J
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "B"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        for I = Int64(v_["2"]):Int64(v_["PT"])
            v_["I-1"] = -1+I
            for J = Int64(v_["2"]):Int64(v_["PY"])
                v_["J-1"] = -1+J
                ename = "C"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(Int64(v_["I-1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "D"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(I)*","*string(Int64(v_["J-1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["PT-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HT"]
            v_["CXI1"] = cos(v_["XI1"])
            v_["ECX"] = v_["CXI1"]*v_["EX"]
            v_["ECX1"] = 1.0+v_["ECX"]
            v_["E12"] = v_["ECX1"]*v_["ECX1"]
            v_["WI"] = v_["ECX1"]*v_["E12"]
            v_["2WI"] = v_["WI"]+v_["WI"]
            v_["XI+1"] = v_["XI1"]+v_["HT"]
            v_["CXI+1"] = cos(v_["XI+1"])
            v_["E+CX0"] = v_["CXI+1"]*v_["EX"]
            v_["E+CX1"] = 1.0+v_["E+CX0"]
            v_["E22"] = v_["E+CX1"]*v_["E+CX1"]
            v_["WI+1"] = v_["E+CX1"]*v_["E22"]
            v_["PM0"] = v_["2WI"]+v_["WI+1"]
            v_["PM1"] = v_["PM0"]/v_["SIX"]
            v_["LA/HY2"] = v_["PM1"]*v_["HT/HY"]
            v_["LA/HT2"] = v_["PM1"]*v_["HY/HT"]
            for J = Int64(v_["1"]):Int64(v_["PY-1"])
                ig = ig_["GR"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["LA/HT2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["LA/HY2"]))
            end
        end
        for I = Int64(v_["2"]):Int64(v_["PT"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["XI1"] = v_["RI-1"]*v_["HT"]
            v_["CXI1"] = cos(v_["XI1"])
            v_["ECX"] = v_["CXI1"]*v_["EX"]
            v_["ECX1"] = 1.0+v_["ECX"]
            v_["E12"] = v_["ECX1"]*v_["ECX1"]
            v_["WI"] = v_["ECX1"]*v_["E12"]
            v_["2WI"] = v_["WI"]+v_["WI"]
            v_["XI-1"] = v_["XI1"]-v_["HT"]
            v_["CXI-1"] = cos(v_["XI-1"])
            v_["E-CX0"] = v_["CXI-1"]*v_["EX"]
            v_["E-CX1"] = 1.0+v_["E-CX0"]
            v_["E32"] = v_["E-CX1"]*v_["E-CX1"]
            v_["WI-1"] = v_["E-CX1"]*v_["E32"]
            v_["PL0"] = v_["2WI"]+v_["WI-1"]
            v_["PL1"] = v_["PL0"]/v_["SIX"]
            v_["MU/HY2"] = v_["PL1"]*v_["HT/HY"]
            v_["MU/HT2"] = v_["PL1"]*v_["HY/HT"]
            for J = Int64(v_["2"]):Int64(v_["PY"])
                ig = ig_["GL"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["MU/HT2"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["D"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["MU/HY2"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(4)            -0.2247400
# LO SOLTN(10)           -0.1789600
# LO SOLTN(23)           -0.1800500
# LO SOLTN(32)           -0.1803000
# LO SOLTN(75)           -0.1805500
# LO SOLTN(100)          -0.1805700
# LO SOLTN(125)          -0.1805800
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-QBR2-AY-V-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

