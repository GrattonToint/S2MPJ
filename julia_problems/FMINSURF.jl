function FMINSURF(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : FMINSURF
#    *********
# 
#    The free boundary minimum surface problem.
# 
#    The problem comes from the discretization of the minimum surface
#    problem on the unit square with "free boundary conditions"
#    one must find the minumum surface over the unit square 
#    (which is clearly 1.0).  Furthermore, the average distance of the surface
#    from zero is also minimized.
# 
#    The Hessian is dense.
# 
#    The unit square is discretized into (p-1)**2 little squares. The
#    heights of the considered surface above the corners of these little
#    squares are the problem variables,  There are p**2 of them.
#    Given these heights, the area above a little square is
#    approximated by the
#      S(i,j) = sqrt( 1 + 0.5(p-1)**2 ( a(i,j) + b(i,j) ) ) / (p-1)**2
#    where
#      a(i,j) = x(i,j) - x(i+1,j+1)
#    and
#      b(i,j) = x(i+1,j) - x(i,j+1)
# 
#    Source: setting the boundary free in 
#    A Griewank and Ph. Toint,
#    "Partitioned variable metric updates for large structured
#    optimization problems",
#    Numerische Mathematik 39:429-448, 1982.
# 
#    SIF input: Ph. Toint, November 1991.
# 
#    classification = "C-OUR2-MY-V-0"
# 
#    P is the number of points in one side of the unit square
# 
#       Alternative values for the SIF file parameters:
# IE P                   4              $-PARAMETER n = 16     original value
# IE P                   7              $-PARAMETER n = 49
# IE P                   8              $-PARAMETER n = 64
# IE P                   11             $-PARAMETER n = 121
# IE P                   31             $-PARAMETER n = 961
# IE P                   32             $-PARAMETER n = 1024
# IE P                   75             $-PARAMETER n = 5625
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "FMINSURF"

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
            v_["P"] = Int64(4);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
# IE P                   100            $-PARAMETER n = 10000
# IE P                   125            $-PARAMETER n = 15625
        v_["H00"] = 1.0
        v_["SLOPEJ"] = 4.0
        v_["SLOPEI"] = 8.0
        v_["TWOP"] = v_["P"]+v_["P"]
        v_["P-1"] = -1+v_["P"]
        v_["PP-1"] = v_["P"]*v_["P-1"]
        v_["RP-1"] = Float64(v_["P-1"])
        v_["INVP-1"] = 1.0/v_["RP-1"]
        v_["RP-1SQ"] = v_["INVP-1"]*v_["INVP-1"]
        v_["SCALE"] = 1.0/v_["RP-1SQ"]
        v_["SQP-1"] = v_["RP-1"]*v_["RP-1"]
        v_["PARAM"] = 0.5*v_["SQP-1"]
        v_["RP"] = Float64(v_["P"])
        v_["P2"] = v_["RP"]*v_["RP"]
        v_["P4"] = v_["P2"]*v_["P2"]
        v_["1"] = 1
        v_["2"] = 2
        v_["STON"] = v_["INVP-1"]*v_["SLOPEI"]
        v_["WTOE"] = v_["INVP-1"]*v_["SLOPEJ"]
        v_["H01"] = v_["H00"]+v_["SLOPEJ"]
        v_["H10"] = v_["H00"]+v_["SLOPEI"]
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
        for I = Int64(v_["1"]):Int64(v_["P-1"])
            for J = Int64(v_["1"]):Int64(v_["P-1"])
                ig,ig_,_ = s2mpj_ii("S"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                arrset(pbm.gscale,ig,Float64(v_["SCALE"]))
            end
        end
        for J = Int64(v_["1"]):Int64(v_["P"])
            for I = Int64(v_["1"]):Int64(v_["P"])
                ig,ig_,_ = s2mpj_ii("AVH",ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        ig,ig_,_ = s2mpj_ii("AVH",ig_)
        arrset(gtype,ig,"<>")
        arrset(pbm.gscale,ig,Float64(v_["P4"]))
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = fill(-1.0,ngrp)
        pbm.gconst[ig_["AVH"]] = Float64(0.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        for J = Int64(v_["1"]):Int64(v_["P"])
            v_["J-1"] = -1+J
            v_["RJ-1"] = Float64(v_["J-1"])
            v_["TH"] = v_["RJ-1"]*v_["WTOE"]
            v_["TL"] = v_["TH"]+v_["H00"]
            v_["TU"] = v_["TH"]+v_["H10"]
            pb.x0[ix_["X"*string(Int64(v_["1"]))*","*string(J)]] = Float64(v_["TL"])
            pb.x0[ix_["X"*string(Int64(v_["P"]))*","*string(J)]] = Float64(v_["TU"])
        end
        for I = Int64(v_["2"]):Int64(v_["P-1"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TV"] = v_["RI-1"]*v_["STON"]
            v_["TR"] = v_["TV"]+v_["H00"]
            v_["TL"] = v_["TV"]+v_["H01"]
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["P"]))]] = Float64(v_["TL"])
            pb.x0[ix_["X"*string(I)*","*string(Int64(v_["1"]))]] = Float64(v_["TR"])
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
        for I = Int64(v_["1"]):Int64(v_["P-1"])
            v_["I+1"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["P-1"])
                v_["J+1"] = 1+J
                ename = "A"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["I+1"]))*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "B"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "X"*string(Int64(v_["I+1"]))*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(I)*","*string(Int64(v_["J+1"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQROOT",igt_)
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["AVH"]
        arrset(pbm.grftype,ig,"gL2")
        for I = Int64(v_["1"]):Int64(v_["P-1"])
            for J = Int64(v_["1"]):Int64(v_["P-1"])
                ig = ig_["S"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gSQROOT")
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["PARAM"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["PARAM"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-MY-V-0"
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQROOT"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SQRAL = sqrt(GVAR_)
        f_= SQRAL
        if nargout>1
            g_ = 0.5e0/SQRAL
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -0.25e0/(SQRAL*GVAR_)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

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
                H_ = 2.0e0
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

