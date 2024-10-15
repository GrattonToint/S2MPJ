function KISSING(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: KISSING NUMBER PROBLEM
#                                                                    
#    Source: This problem is associated to the family of Hard-Spheres 
#    problem. It belongs to the family of sphere packing problems, a 
#    class of challenging problems dating from the beginning of the 
#    17th century which is related to practical problems in Chemistry, 
#    Biology and Physics. It consists on maximizing the minimum pairwise 
#    distance between NP points on a sphere in \R^{MDIM}. 
#    This problem may be reduced to a nonconvex nonlinear optimization 
#    problem with a potentially large number of (nonoptimal) points 
#    satisfying optimality conditions. We have, thus, a class of problems 
#    indexed by the parameters MDIM and NP, that provides a suitable 
#    set of test problems for evaluating nonlinear programming codes.
#    After some algebric manipulations, we can formulate this problem as
#                             Minimize z
#                             subject to
#        
#       z \geq <x_i, x_j> for all different pair of indices i, j
#       
#                             ||x_i||^2 = 1    for all i = 1,...,NP
#      The goal is to find an objective value less than 0.5 (This means
#      that the NP points stored belong to the sphere and every distance
#      between two of them is greater than 1.0).
#      Obs: the starting point is aleatorally chosen although each 
#      variable belongs to [-1.,1.].
#      References:
#      [1] "Validation of an Augmented Lagrangian algorithm with a 
#           Gauss-Newton Hessian approximation using a set of 
#           Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello 
#           and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP, 
#           Campinas, 1998.
#      [2] "Inexact-Restoration Algorithm for Constrained Optimization",
#           J. M. Martinez and E. A. Pilotta, Tech. Report, IMECC-UNICAMP, 
#           Campinas, 1998.
#      [3]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
#            N. J. C. Sloane, Springer-Verlag, NY, 1988.
#      SIF input: September 29, 1998
# 		 Jose Mario Martinez
#                 Elvio Angel Pilotta
# 
#    classification = "C-LQR2-RN-V-V"
# 
# **********************************************************************
# 
#    Number of points: NP >= 12
# 
#       Alternative values for the SIF file parameters:
# IE NP                   12            $-PARAMETER
# IE NP                   13            $-PARAMETER
# IE NP                   14            $-PARAMETER
# IE NP                   15            $-PARAMETER
# IE NP                   22            $-PARAMETER
# IE NP                   23            $-PARAMETER
# IE NP                   24            $-PARAMETER
# IE NP                   25            $-PARAMETER
# IE NP                   26            $-PARAMETER
# IE NP                   27            $-PARAMETER
# IE NP	                 37            $-PARAMETER
# IE NP                   38            $-PARAMETER
# IE NP                   39            $-PARAMETER
# IE NP                   40            $-PARAMETER
# IE NP                   41            $-PARAMETER
# IE NP                   42            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "KISSING"

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
            v_["NP"] = Int64(25);  #  SIF file default value
        else
            v_["NP"] = Int64(args[1]);
        end
# IE MDIM                 3             $-PARAMETER
        if nargin<2
            v_["MDIM"] = Int64(3);  #  SIF file default value
        else
            v_["MDIM"] = Int64(args[2]);
        end
# IE MDIM                 4             $-PARAMETER
# IE MDIM                 5             $-PARAMETER
        v_["N-"] = -1+v_["NP"]
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NP"])
            for J = Int64(v_["1"]):Int64(v_["MDIM"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
            end
        end
        iv,ix_,_ = s2mpj_ii("Z",ix_)
        arrset(pb.xnames,iv,"Z")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Z"]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["N-"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP"])
                ig,ig_,_ = s2mpj_ii("IC"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<=")
                arrset(pb.cnames,ig,"IC"*string(I)*","*string(J))
                iv = ix_["Z"]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NP"])
            ig,ig_,_ = s2mpj_ii("EC"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EC"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["NP"])
            pbm.gconst[ig_["EC"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NP"])
            for J = Int64(v_["1"]):Int64(v_["MDIM"])
                pb.xlower[ix_["X"*string(I)*","*string(J)]] = -Inf
                pb.xupper[ix_["X"*string(I)*","*string(J)]] = +Inf
            end
        end
        pb.xlower[ix_["Z"]] = -Inf
        pb.xupper[ix_["Z"]] = +Inf
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X1,1"]] = Float64(-0.10890604)
        pb.x0[ix_["X1,2"]] = Float64(0.85395078)
        pb.x0[ix_["X1,3"]] = Float64(-0.45461680)
        pb.x0[ix_["X2,1"]] = Float64(0.49883922)
        pb.x0[ix_["X2,2"]] = Float64(-0.18439316)
        pb.x0[ix_["X2,3"]] = Float64(-0.04798594)
        pb.x0[ix_["X3,1"]] = Float64(0.28262888)
        pb.x0[ix_["X3,2"]] = Float64(-0.48054070)
        pb.x0[ix_["X3,3"]] = Float64(0.46715332)
        pb.x0[ix_["X4,1"]] = Float64(-0.00580106)
        pb.x0[ix_["X4,2"]] = Float64(-0.49987584)
        pb.x0[ix_["X4,3"]] = Float64(-0.44130302)
        pb.x0[ix_["X5,1"]] = Float64(0.81712540)
        pb.x0[ix_["X5,2"]] = Float64(-0.36874258)
        pb.x0[ix_["X5,3"]] = Float64(-0.68321896)
        pb.x0[ix_["X6,1"]] = Float64(0.29642426)
        pb.x0[ix_["X6,2"]] = Float64(0.82315508)
        pb.x0[ix_["X6,3"]] = Float64(0.35938150)
        pb.x0[ix_["X7,1"]] = Float64(0.09215152)
        pb.x0[ix_["X7,2"]] = Float64(-0.53564686)
        pb.x0[ix_["X7,3"]] = Float64(0.00191436)
        pb.x0[ix_["X8,1"]] = Float64(0.11700318)
        pb.x0[ix_["X8,2"]] = Float64(0.96722760)
        pb.x0[ix_["X8,3"]] = Float64(-0.14916438)
        pb.x0[ix_["X9,1"]] = Float64(0.01791524)
        pb.x0[ix_["X9,2"]] = Float64(0.17759446)
        pb.x0[ix_["X9,3"]] = Float64(-0.61875872)
        pb.x0[ix_["X10,1"]] = Float64(-0.63833630)
        pb.x0[ix_["X10,2"]] = Float64(0.80830972)
        pb.x0[ix_["X10,3"]] = Float64(0.45846734)
        pb.x0[ix_["X11,1"]] = Float64(0.28446456)
        pb.x0[ix_["X11,2"]] = Float64(0.45686938)
        pb.x0[ix_["X11,3"]] = Float64(0.16368980)
        pb.x0[ix_["X12,1"]] = Float64(0.76557382)
        pb.x0[ix_["X12,2"]] = Float64(0.16700944)
        pb.x0[ix_["X12,3"]] = Float64(-0.31647534)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eQUA", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N-"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP"])
                for K = Int64(v_["1"]):Int64(v_["MDIM"])
                    ename = "A"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"ePROD")
                    arrset(ielftype,ie,iet_["ePROD"])
                    vname = "X"*string(I)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "X"*string(J)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NP"])
            for K = Int64(v_["1"]):Int64(v_["MDIM"])
                ename = "B"*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eQUA")
                arrset(ielftype,ie,iet_["eQUA"])
                vname = "X"*string(I)*","*string(K)
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
        for I = Int64(v_["1"]):Int64(v_["N-"])
            v_["I+"] = 1+I
            for J = Int64(v_["I+"]):Int64(v_["NP"])
                for K = Int64(v_["1"]):Int64(v_["MDIM"])
                    ig = ig_["IC"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["A"*string(I)*","*string(J)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["NP"])
            for K = Int64(v_["1"]):Int64(v_["MDIM"])
                ig = ig_["EC"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             4.47214D-01
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-RN-V-V"
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

    elseif action == "ePROD"

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

    elseif action == "eQUA"

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

