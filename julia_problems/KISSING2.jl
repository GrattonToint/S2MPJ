function KISSING2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem: A second formulation of the KISSING NUMBER PROBLEM
#                                                                    
#    Source: This problem is associated to the family of Hard-Spheres 
#    problem. It belongs to the family of sphere packing problems, a 
#    class of challenging problems dating from the beginning of the 
#    17th century which is related to practical problems in Chemistry, 
#    Biology and Physics. Given a fixed unit sphere at the origin in R^n, 
#    the problem consists of arranging a further m unit spheres so that 
#    sum of the distances to these spheres is as small as possible.
#    This problem may be reduced to a nonconvex nonlinear optimization 
#    problem with a potentially large number of (nonoptimal) points 
#    satisfying optimality conditions. We have, thus, a class of problems 
#    indexed by the parameters m and n, that provides a suitable 
#    set of test problems for evaluating nonlinear programming codes.
#    After some algebric manipulations, we can formulate this problem as
#               m
#     Minimize sum <p_i,p_i> - m n
#              i=1
#     subject to
#        
#      <p_i - p_j, p_i - p_j> >= 4 for all different pair of indices i, j
#      and  
#      <p_i, p_i> >= 4 for all indices i
#   
#      as well as n(n-1)/2 normalisation constraints fixing components.
#      The goal is to find an objective value equal to 0.
#      [1]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
#            N. J. C. Sloane, Springer-Verlag, NY, 1988.
#    SIF input: Nick Gould, September 2000
# 
#    classification = "C-QQR2-RN-V-V"
# 
# **********************************************************************
# 
#    Number of points: m
# 
#       Alternative values for the SIF file parameters:
# IE m                   24             $-PARAMETER  number of points
# IE m                   25             $-PARAMETER  number of points
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "KISSING2"

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
            v_["m"] = Int64(25);  #  SIF file default value
        else
            v_["m"] = Int64(args[1]);
        end
# IE m                   100            $-PARAMETER  number of points
# IE n                    4             $-PARAMETER  dimension of sphere
        if nargin<2
            v_["n"] = Int64(4);  #  SIF file default value
        else
            v_["n"] = Int64(args[2]);
        end
# IE n                    8             $-PARAMETER  dimension of sphere
        v_["1"] = 1
        v_["2"] = 2
        v_["n-1"] = v_["n"]-v_["1"]
        v_["rm"] = Float64(v_["m"])
        v_["rn"] = Float64(v_["n"])
        v_["RM+N"] = v_["rm"]+v_["rn"]
        v_["mn"] = v_["rm"]*v_["rn"]
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["PI/m"] = v_["PI"]/v_["rm"]
        v_["2PI/m"] = 2.0*v_["PI/m"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["m"])
            for J = Int64(v_["1"]):Int64(v_["n"])
                iv,ix_,_ = s2mpj_ii("P"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"P"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["m"])
            for J = Int64(v_["1"]):Int64(v_["m"])
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
        pbm.gconst[ig_["OBJ"]] = Float64(v_["mn"])
        for I = Int64(v_["1"]):Int64(v_["m"])
            for J = Int64(v_["1"]):Int64(v_["m"])
                pbm.gconst[ig_["C"*string(I)*","*string(J)]] = Float64(4.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["m"])
            for J = Int64(v_["1"]):Int64(v_["n"])
                pb.xlower[ix_["P"*string(I)*","*string(J)]] = -Inf
                pb.xupper[ix_["P"*string(I)*","*string(J)]] = +Inf
            end
        end
        for I = Int64(v_["2"]):Int64(v_["n"])
            for J = Int64(I):Int64(v_["n"])
                pb.xlower[ix_["P"*string(I)*","*string(J)]] = 0.0
                pb.xupper[ix_["P"*string(I)*","*string(J)]] = 0.0
            end
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["m"])
            v_["RI"] = Float64(I)
            v_["2PIi/m"] = v_["2PI/m"]*v_["RI"]
            v_["cos"] = cos(v_["2PIi/m"])
            v_["sin"] = sin(v_["2PIi/m"])
            v_["cos"] = v_["cos"]+v_["cos"]
            v_["sin"] = v_["sin"]+v_["sin"]
            if haskey(ix_,"P"*string(I)*","*string(Int64(v_["1"])))
                pb.x0[ix_["P"*string(I)*","*string(Int64(v_["1"]))]] = Float64(v_["cos"])
            else
                pb.y0[findfirst(x->x==ig_["P"*string(I)*","*string(Int64(v_["1"]))],pbm.congrps)] = Float64(v_["cos"])
            end
            for J = Int64(v_["2"]):Int64(v_["n-1"])
                if haskey(ix_,"P"*string(I)*","*string(J))
                    pb.x0[ix_["P"*string(I)*","*string(J)]] = Float64(v_["sin"])
                else
                    pb.y0[findfirst(x->x==ig_["P"*string(I)*","*string(J)],pbm.congrps)]  = (
                          Float64(v_["sin"]))
                end
            end
        end
        if haskey(ix_,"P"*string(Int64(v_["m"]))*","*string(Int64(v_["n"])))
            pb.x0[ix_["P"*string(Int64(v_["m"]))*","*string(Int64(v_["n"]))]]  = (
                  Float64(v_["cos"]))
        else
            pb.y0[findfirst(x->x==ig_["P"*string(Int64(v_["m"]))*","*string(Int64(v_["n"]))],pbm.congrps)] = Float64(v_["cos"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD1", iet_)
        loaset(elftv,it,1,"P")
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"Q")
        loaset(elftv,it,2,"R")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["m"])
            v_["I-"] = -1+I
            v_["I+"] = 1+I
            for J = Int64(v_["1"]):Int64(v_["I-"])
                for K = Int64(v_["1"]):Int64(v_["n"])
                    ename = "E"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"ePROD2")
                    arrset(ielftype,ie,iet_["ePROD2"])
                    vname = "P"*string(I)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="Q",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "P"*string(J)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="R",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
            for K = Int64(v_["1"]):Int64(v_["n"])
                ename = "E"*string(I)*","*string(I)*","*string(K)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"ePROD1")
                arrset(ielftype,ie,iet_["ePROD1"])
                vname = "P"*string(I)*","*string(K)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="P",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for J = Int64(v_["I+"]):Int64(v_["m"])
                for K = Int64(v_["1"]):Int64(v_["n"])
                    ename = "E"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,_  = s2mpj_ii(ename,ie_)
                    arrset(pbm.elftype,ie,"ePROD2")
                    arrset(ielftype,ie,iet_["ePROD2"])
                    vname = "P"*string(I)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="Q",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "P"*string(J)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="R",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["m"])
            for K = Int64(v_["1"]):Int64(v_["n"])
                ig = ig_["OBJ"]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(I)*","*string(K)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
            for J = Int64(v_["1"]):Int64(v_["m"])
                for K = Int64(v_["1"]):Int64(v_["n"])
                    ig = ig_["C"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION             0.00000D+00   $ n=4, m = 24
# XL SOLUTION             6.48030D+00   $ n=4, m = 25 one of many local solutions
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QQR2-RN-V-V"
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

    elseif action == "ePROD1"

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

    elseif action == "ePROD2"

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

