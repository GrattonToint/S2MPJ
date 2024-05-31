function POWELLBC(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POWELLBC
#    --------
# 
#    A bound-constrained optimization problem to 
#    separate points within a square in the plane
# 
#    Given points p_j = ( x_2j-1 , x_2j ), i = 1, ..., p
# 
#    minimize sum_k=2^p sum_j=1^k-1 1 / || p_j - p_k ||_2
# 
#    subject to 0 <= x_i <= 1, i = 1, ..., 2p = n
# 
#    Source: 
#    M. J. D. Powell
#    Private communication (Optbridge, 2006)
# 
#    SIF input: Nick Gould, Aug 2006.
# 
#    classification = "OBR2-AN-V-0"
# 
#    Number of points
# 
#       Alternative values for the SIF file parameters:
# IE P                   2              $-PARAMETER
# IE P                   5              $-PARAMETER
# IE P                   10             $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "POWELLBC"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "POWELLBC"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["P"] = Int64(12);  #  SIF file default value
        else
            v_["P"] = Int64(args[1]);
        end
# IE P                   100            $-PARAMETER
# IE P                   500            $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["N"] = 2*v_["P"]
        v_["RN"] = Float64(v_["N"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.0,pb.n)
        pb.xupper = fill(1.0,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["T"] = v_["RI"]/v_["RN"]
            v_["T"] = v_["T"]*v_["T"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["T"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eINVNRM", iet_)
        loaset(elftv,it,1,"XJ")
        loaset(elftv,it,2,"YJ")
        loaset(elftv,it,3,"XK")
        loaset(elftv,it,4,"YK")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for K = Int64(v_["2"]):Int64(v_["P"])
            v_["K-1"] = K-v_["1"]
            v_["2K"] = v_["2"]*K
            v_["2K-1"] = v_["2K"]-v_["1"]
            for J = Int64(v_["1"]):Int64(v_["K-1"])
                v_["2J"] = v_["2"]*J
                v_["2J-1"] = v_["2J"]-v_["1"]
                ename = "E"*string(K)*","*string(J)
                ie,ie_,newelt = s2x_ii(ename,ie_)
                if newelt > 0
                    arrset(pbm.elftype,ie,"eINVNRM")
                    arrset(ielftype,ie,iet_["eINVNRM"])
                end
                vname = "X"*string(Int64(v_["2J-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,0.0,1.0,nothing)
                posev = findfirst(x->x=="XJ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["2K-1"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,0.0,1.0,nothing)
                posev = findfirst(x->x=="XK",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["2J"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,0.0,1.0,nothing)
                posev = findfirst(x->x=="YJ",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(Int64(v_["2K"]))
                iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,0.0,1.0,nothing)
                posev = findfirst(x->x=="YK",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for K = Int64(v_["2"]):Int64(v_["P"])
            v_["K-1"] = K-v_["1"]
            for J = Int64(v_["1"]):Int64(v_["K-1"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(K)*","*string(J)])
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "OBR2-AN-V-0"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eINVNRM"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        NORM = 1.0/sqrt(IV_[1]*IV_[1]+IV_[2]*IV_[2])
        f_   = NORM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -IV_[1]*NORM^3
            g_[2] = -IV_[2]*NORM^3
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = (3.0*IV_[1]*IV_[1]*NORM^2-1.0)*NORM^3
                H_[1,2] = 3.0*IV_[1]*IV_[2]*NORM^5
                H_[2,1] = H_[1,2]
                H_[2,2] = (3.0*IV_[2]*IV_[2]*NORM^2-1.0)*NORM^3
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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

