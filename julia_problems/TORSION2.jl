function TORSION2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TORSION2
#    *********
# 
#    The quadratic elastic torsion problem
# 
#    The problem comes from the obstacle problem on a square.
# 
#    The square is discretized into (px-1)(py-1) little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are px**2 of them.
# 
#    The dimension of the problem is specified by Q, which is half the
#    number discretization points along one of the coordinate
#    direction.
#    Since the number of variables is P**2, it is given by 4Q**2
# 
#    Source: problem (c=5, starting point Z = origin) in
#    J. More' and G. Toraldo,
#    "On the Solution of Large Quadratic-Programming Problems with Bound
#    Constraints", 
#    SIAM J. on Optimization, vol 1(1), pp. 93-113, 1991.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-QBR2-MY-V-0"
# 
#       Alternative values for the SIF file parameters:
# IE Q                   2              $-PARAMETER n= 16
# IE Q                   5              $-PARAMETER n= 100     original value
# IE Q                   11             $-PARAMETER n= 484
# IE Q                   16             $-PARAMETER n= 1024
# IE Q                   37             $-PARAMETER n= 5476
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TORSION2"

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
            v_["Q"] = Int64(2);  #  SIF file default value
        else
            v_["Q"] = Int64(args[1]);
        end
# IE Q                   50             $-PARAMETER n= 10000
# IE Q                   61             $-PARAMETER n= 14884
        if nargin<2
            v_["C"] = Float64(5.0);  #  SIF file default value
        else
            v_["C"] = Float64(args[2]);
        end
        v_["Q+1"] = 1+v_["Q"]
        v_["P"] = v_["Q"]+v_["Q"]
        v_["P-1"] = -1+v_["P"]
        v_["1/H"] = Float64(v_["P-1"])
        v_["H"] = 1.0/v_["1/H"]
        v_["H2"] = v_["H"]*v_["H"]
        v_["C0"] = v_["H2"]*v_["C"]
        v_["LC"] = -1.0*v_["C0"]
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            for J = Int64(v_["2"]):Int64(v_["P-1"])
                ig,ig_,_ = s2mpj_ii("G"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(v_["LC"])
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
        for J = Int64(v_["1"]):Int64(v_["P"])
            pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = 0.0
            pb.xlower[ix_["X"*string(Int64(v_["P"]))*","*string(J)]] = 0.0
            pb.xupper[ix_["X"*string(Int64(v_["P"]))*","*string(J)]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["P"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["P"]))]] = 0.0
            pb.xlower[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
            pb.xupper[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = 0.0
        end
        for I = Int64(v_["2"]):Int64(v_["Q"])
            for J = Int64(v_["2"]):Int64(I)
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["UPPL"] = v_["RJ-1"]*v_["H"]
                v_["LOWL"] = -1.0*v_["UPPL"]
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWL"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPL"]
            end
            v_["MI"] = -1*I
            v_["P-I"] = v_["P"]+v_["MI"]
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["UPPM"] = v_["RI-1"]*v_["H"]
            v_["LOWM"] = -1.0*v_["UPPM"]
            v_["P-I+1"] = 1+v_["P-I"]
            for J = Int64(I):Int64(v_["P-I+1"])
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWM"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPM"]
            end
            for J = Int64(v_["P-I+1"]):Int64(v_["P-1"])
                v_["MJ"] = -1*J
                v_["P-J"] = v_["P"]+v_["MJ"]
                v_["RP-J"] = Float64(v_["P-J"])
                v_["UPPR"] = v_["RP-J"]*v_["H"]
                v_["LOWR"] = -1.0*v_["UPPR"]
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWR"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPR"]
            end
        end
        for I = Int64(v_["Q+1"]):Int64(v_["P-1"])
            v_["MI"] = -1*I
            v_["P-I"] = v_["P"]+v_["MI"]
            v_["P-I+1"] = 1+v_["P-I"]
            for J = Int64(v_["2"]):Int64(v_["P-I+1"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["UPPL"] = v_["RJ-1"]*v_["H"]
                v_["LOWL"] = -1.0*v_["UPPL"]
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWL"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPL"]
            end
            v_["RP-I"] = Float64(v_["P-I"])
            v_["UPPM"] = v_["RP-I"]*v_["H"]
            v_["LOWM"] = -1.0*v_["UPPM"]
            for J = Int64(v_["P-I+1"]):Int64(I)
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWM"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPM"]
            end
            for J = Int64(I):Int64(v_["P-1"])
                v_["MJ"] = -1*J
                v_["P-J"] = v_["P"]+v_["MJ"]
                v_["RP-J"] = Float64(v_["P-J"])
                v_["UPPR"] = v_["RP-J"]*v_["H"]
                v_["LOWR"] = -1.0*v_["UPPR"]
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = v_["LOWR"]
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = v_["UPPR"]
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
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            for J = Int64(v_["2"]):Int64(v_["P-1"])
                v_["J-1"] = -1+J
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
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            for J = Int64(v_["2"]):Int64(v_["P-1"])
                ig = ig_["G"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(0.25))
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(0.25))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(0.25))
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["D"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(0.25))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(2)            -5.1851852D-1
# LO SOLTN(5)            -4.9234185D-1
# LO SOLTN(11)           -4.5608771D-1
# LO SOLTN(16)           ???
# LO SOLTN(37)           ???
# LO SOLTN(50)           ???
# LO SOLTN(61)           ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-QBR2-MY-V-0"
        pb.x0          = zeros(Float64,pb.n)
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

