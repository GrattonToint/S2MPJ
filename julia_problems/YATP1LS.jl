function YATP1LS(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : YATP1LS
#    *********
# 
#    Yet another test problem involving double pseudo-stochastic constraints
#    on a square matrix. If the matrix dimension is N, the number of
#    variables is equal to  N**2 + 2* N. The equations are
#    x_{ij}^3 - A x_{ij}^2 - ( y_i + z_i ) ( x_{ij}cos(x_{ij} - sin(x_{ij}) ) = 0
#                                                     (i,j = 1, ..., N )
#    \sum_i^N sin(x_{ij}) / x_{ij} = 1                  (j = 1,..., N)
#    \sum_j^N sin(x_{ij}) / x_{ij} = 1                  (i = 1,..., N)
#    The problem is non convex.
# 
#   Least-squares version
# 
#    Source:
#    a late evening idea by Ph. Toint
# 
#    SIF input: Ph. Toint, June 2003.
# 
#    classification = "C-CSUR2-AN-V-V"
# 
#    The dimension of the matrix
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n = 120
# IE N                   50             $-PARAMETER n = 2600
# IE N                   100            $-PARAMETER n = 10200
# IE N                   200            $-PARAMETER n = 40400
# IE N                   350            $-PARAMETER n = 123200
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "YATP1LS"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling YATP1LS.")
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
            v_["N"] = Int64(5);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
        v_["A"] = 10.0
        v_["1"] = 1
        v_["-A"] = -1.0*v_["A"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2mpj_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("E"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("ER"*string(I),ig_)
            arrset(gtype,ig,"<>")
            ig,ig_,_ = s2mpj_ii("EC"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["ER"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["EC"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                pb.x0[ix_["X"*string(I)*","*string(J)]] = Float64(6.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eCB", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eLXC", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eLS", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"Z")
        it,iet_,_ = s2mpj_ii( "eRAT", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ename = "CB"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eCB")
                arrset(ielftype,ie,iet_["eCB"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "LS"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DC"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eLXC")
                arrset(ielftype,ie,iet_["eLXC"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Z"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DS"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eLS")
                arrset(ielftype,ie,iet_["eLS"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Z"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "SX"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eRAT")
                arrset(ielftype,ie,iet_["eRAT"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["E"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["CB"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(1.0))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["LS"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["-A"]))
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["DC"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["DS"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(1.0))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["N"])
                ig = ig_["ER"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["SX"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
                ig = ig_["EC"*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["SX"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-CSUR2-AN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm


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

    elseif action == "eCB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*EV_[1]*EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 6.0*EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eLXC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        C = cos(IV_[1])
        S = sin(IV_[1])
        f_   = IV_[2]*IV_[1]*C
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = IV_[1]*C
            g_[1] = IV_[2]*(C-IV_[1]*S)
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = C-IV_[1]*S
                H_[1,2] = H_[2,1]
                H_[1,1] = -IV_[2]*(S+S+IV_[1]*C)
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

    elseif action == "eLS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        C = cos(IV_[1])
        S = sin(IV_[1])
        f_   = IV_[2]*S
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = S
            g_[1] = IV_[2]*C
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = C
                H_[1,2] = H_[2,1]
                H_[1,1] = -IV_[2]*S
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

    elseif action == "eRAT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        C = cos(EV_[1])
        S = sin(EV_[1])
        f_   = S/EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (C-S/EV_[1])/EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -S/EV_[1]-(C+C)/EV_[1]^2+(S+S)/EV_[1]^3
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

