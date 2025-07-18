function QR3D(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : QR3D
#    *********
# 
#    Find the QR factorization of a tridiagonal matrix A.
#    The problem is formulated as a system of quadratic equations
#    whose unknowns are the elements of the orthogonal matrix Q and of
#    the upper triangular matrix R.  In this version of the problem,
#    the banded structure of R is not imposed as a constraint. See problem
#    QR3DBD for the case where this structure is explicitly used.
# 
#    The problem is non-convex.
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, Nov 1993
# 
#    classification = "C-CNQR2-AN-V-V"
# 
#    Define the matrix order M  ( M >= 3 ).
#    There are M * ( 3M + 1) / 2 variables and equations.
# 
#       Alternative values for the SIF file parameters:
# IE M                   5              $-PARAMETER  n =  40
# IE M                   10             $-PARAMETER  n = 155  original value
# IE M                   20             $-PARAMETER  n = 610
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "QR3D"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling QR3D.")
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
            v_["M"] = Int64(5);  #  SIF file default value
        else
            v_["M"] = Int64(args[1]);
        end
        v_["1"] = 1
        v_["2"] = 2
        v_["M-1"] = -1+v_["M"]
        v_["RM"] = Float64(v_["M"])
        v_["2/M"] = 2.0/v_["RM"]
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))] = v_["2/M"]
        v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))] = 0.0
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["1-I"] = -1*v_["I-1"]
            v_["R1-I"] = Float64(v_["1-I"])
            v_["1-I/M"] = v_["R1-I"]/v_["RM"]
            v_["2I"] = 2*I
            v_["R2I"] = Float64(v_["2I"])
            v_["2I/M"] = v_["R2I"]/v_["RM"]
            v_["A"*string(I)*","*string(Int64(v_["I-1"]))] = v_["1-I/M"]
            v_["A"*string(I)*","*string(I)] = v_["2I/M"]
            v_["A"*string(I)*","*string(Int64(v_["I+1"]))] = v_["1-I/M"]
        end
        v_["RM-1"] = Float64(v_["M-1"])
        v_["1-M"] = -1.0*v_["RM-1"]
        v_["1-M/M"] = v_["1-M"]/v_["RM"]
        v_["2M"] = 2.0*v_["RM"]
        v_["A"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))] = v_["1-M/M"]
        v_["A"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))] = v_["2M"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                iv,ix_,_ = s2mpj_ii("Q"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Q"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(I):Int64(v_["M"])
                iv,ix_,_ = s2mpj_ii("R"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"R"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(I):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("O"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"O"*string(I)*","*string(J))
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                ig,ig_,_ = s2mpj_ii("F"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"F"*string(I)*","*string(J))
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["O"*string(I)*","*string(I)]] = Float64(1.0)
        end
        pbm.gconst[ig_["F"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]]  = (
              Float64(v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["1"]))]))
        pbm.gconst[ig_["F"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]]  = (
              Float64(v_["A"*string(Int64(v_["1"]))*","*string(Int64(v_["2"]))]))
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            pbm.gconst[ig_["F"*string(I)*","*string(Int64(v_["I-1"]))]]  = (
                  Float64(v_["A"*string(I)*","*string(Int64(v_["I-1"]))]))
            pbm.gconst[ig_["F"*string(I)*","*string(I)]]  = (
                  Float64(v_["A"*string(I)*","*string(I)]))
            pbm.gconst[ig_["F"*string(I)*","*string(Int64(v_["I+1"]))]]  = (
                  Float64(v_["A"*string(I)*","*string(Int64(v_["I+1"]))]))
        end
        pbm.gconst[ig_["F"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]]  = (
              Float64(v_["A"*string(Int64(v_["M"]))*","*string(Int64(v_["M-1"]))]))
        pbm.gconst[ig_["F"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]]  = (
              Float64(v_["A"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pb.xlower[ix_["R"*string(I)*","*string(I)]] = 0.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pb.x0[ix_["Q"*string(I)*","*string(I)]] = Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["M-1"])
            v_["I+1"] = 1+I
            pb.x0[ix_["R"*string(I)*","*string(I)]]  = (
                  Float64(v_["A"*string(I)*","*string(I)]))
            pb.x0[ix_["R"*string(I)*","*string(Int64(v_["I+1"]))]]  = (
                  Float64(v_["A"*string(I)*","*string(Int64(v_["I+1"]))]))
        end
        pb.x0[ix_["R"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]]  = (
              Float64(v_["A"*string(Int64(v_["M"]))*","*string(Int64(v_["M"]))]))
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(I):Int64(v_["M"])
                for K = Int64(v_["1"]):Int64(v_["M"])
                    ename = "C"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,newelt = s2mpj_ii(ename,ie_)
                    if newelt > 0
                        arrset(pbm.elftype,ie,"en2PR")
                        arrset(ielftype,ie,iet_["en2PR"])
                    end
                    vname = "Q"*string(I)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "Q"*string(J)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                for K = Int64(v_["1"]):Int64(J)
                    ename = "B"*string(I)*","*string(J)*","*string(K)
                    ie,ie_,newelt = s2mpj_ii(ename,ie_)
                    if newelt > 0
                        arrset(pbm.elftype,ie,"en2PR")
                        arrset(ielftype,ie,iet_["en2PR"])
                    end
                    vname = "Q"*string(I)*","*string(K)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                    vname = "R"*string(K)*","*string(J)
                    iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                    posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
                    loaset(pbm.elvar,ie,posev,iv)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(I):Int64(v_["M"])
                for K = Int64(v_["1"]):Int64(v_["M"])
                    ig = ig_["O"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["C"*string(I)*","*string(J)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            for J = Int64(v_["1"]):Int64(v_["M"])
                for K = Int64(v_["1"]):Int64(J)
                    ig = ig_["F"*string(I)*","*string(J)]
                    posel = length(pbm.grelt[ig])+1
                    loaset(pbm.grelt,ig,posel,ie_["B"*string(I)*","*string(J)*","*string(K)])
                    arrset(nlc,length(nlc)+1,ig)
                    loaset(pbm.grelw,ig,posel,1.)
                end
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNQR2-AN-V-V"
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

