function LUKVLE6(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLE6
#    *********
# 
#    Source: Problem 5.6, Generalized Broyden banded function with 
#    exponential constraints, due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-OOR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#       Alternative values for the SIF file parameters:
# IE N                   99             $-PARAMETER
# IE N                   999            $-PARAMETER
# IE N                   9999           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKVLE6"

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
            v_["N"] = Int64(9);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   99999          $-PARAMETER
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["5"] = 5
        v_["6"] = 6
        v_["N/2"] = trunc(Int,(v_["N"]/v_["2"]))
        v_["N+1"] = 1+v_["N"]
        v_["N-4"] = -4+v_["N"]
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
            v_["I+1"] = 1+I
            v_["I-5"] = -5+I
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(2.0)
            v_["A"] = v_["I-5"]
            v_["B"] = v_["1"]
            v_["A"] = Float64(v_["A"])
            v_["ABSA"] = abs(v_["A"])
            v_["ABSA"] = trunc(Int,v_["ABSA"])
            v_["B"] = Float64(v_["B"])
            v_["ABSB"] = abs(v_["B"])
            v_["ABSB"] = trunc(Int,v_["ABSB"])
            v_["ABSA+B"] = v_["ABSA"]+v_["ABSB"]
            v_["A"] = v_["A"]+v_["ABSA+B"]
            v_["B"] = v_["B"]+v_["ABSA+B"]
            v_["A/B"] = trunc(Int,(v_["A"]/v_["B"]))
            v_["B/A"] = trunc(Int,(v_["B"]/v_["A"]))
            v_["SUM"] = v_["A/B"]+v_["B/A"]
            v_["A"] = v_["A"]*v_["A/B"]
            v_["B"] = v_["B"]*v_["B/A"]
            v_["MAXA,B"] = v_["A"]+v_["B"]
            v_["MAXA,B"] = trunc(Int,(v_["MAXA,B"]/v_["SUM"]))
            v_["MAXA,B"] = v_["MAXA,B"]-v_["ABSA+B"]
            v_["MAXI-5,1"] = v_["MAXA,B"]
            v_["A"] = v_["I+1"]
            v_["B"] = v_["N"]
            v_["A"] = -1*v_["A"]
            v_["B"] = -1*v_["B"]
            v_["A"] = Float64(v_["A"])
            v_["ABSA"] = abs(v_["A"])
            v_["ABSA"] = trunc(Int,v_["ABSA"])
            v_["B"] = Float64(v_["B"])
            v_["ABSB"] = abs(v_["B"])
            v_["ABSB"] = trunc(Int,v_["ABSB"])
            v_["ABSA+B"] = v_["ABSA"]+v_["ABSB"]
            v_["A"] = v_["A"]+v_["ABSA+B"]
            v_["B"] = v_["B"]+v_["ABSA+B"]
            v_["A/B"] = trunc(Int,(v_["A"]/v_["B"]))
            v_["B/A"] = trunc(Int,(v_["B"]/v_["A"]))
            v_["SUM"] = v_["A/B"]+v_["B/A"]
            v_["A"] = v_["A"]*v_["A/B"]
            v_["B"] = v_["B"]*v_["B/A"]
            v_["MAXA,B"] = v_["A"]+v_["B"]
            v_["MAXA,B"] = trunc(Int,(v_["MAXA,B"]/v_["SUM"]))
            v_["MINA,B"] = v_["ABSA+B"]-v_["MAXA,B"]
            v_["MINI+1,N"] = v_["MINA,B"]
            for J = Int64(v_["MAXI-5,1"]):Int64(v_["MINI+1,N"])
                ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
            end
        end
        for K = Int64(v_["1"]):Int64(v_["N/2"])
            v_["2K"] = 2*K
            ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(K))
            iv = ix_["X"*string(Int64(v_["2K"]))]
            pbm.A[ig,iv] += Float64(4.0)
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["OBJ"*string(I)]] = Float64(-1.0)
        end
        for K = Int64(v_["1"]):Int64(v_["N/2"])
            pbm.gconst[ig_["C"*string(K)]] = Float64(3.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(3.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eCUBE", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eXEXP", iet_)
        loaset(elftv,it,1,"VM")
        loaset(elftv,it,2,"VP")
        loaset(elftv,it,3,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "S"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCUBE")
            arrset(ielftype,ie,iet_["eCUBE"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for K = Int64(v_["1"]):Int64(v_["N/2"])
            v_["2K"] = 2*K
            v_["2K-1"] = -1+v_["2K"]
            v_["2K+1"] = 1+v_["2K"]
            ename = "P"*string(K)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXEXP")
            arrset(ielftype,ie,iet_["eXEXP"])
            vname = "X"*string(Int64(v_["2K-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VM",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2K"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["2K+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VP",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL7d3",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I+1"] = 1+I
            v_["I-5"] = -5+I
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gL7d3")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(5.0))
            v_["A"] = v_["I-5"]
            v_["B"] = v_["1"]
            v_["A"] = Float64(v_["A"])
            v_["ABSA"] = abs(v_["A"])
            v_["ABSA"] = trunc(Int,v_["ABSA"])
            v_["B"] = Float64(v_["B"])
            v_["ABSB"] = abs(v_["B"])
            v_["ABSB"] = trunc(Int,v_["ABSB"])
            v_["ABSA+B"] = v_["ABSA"]+v_["ABSB"]
            v_["A"] = v_["A"]+v_["ABSA+B"]
            v_["B"] = v_["B"]+v_["ABSA+B"]
            v_["A/B"] = trunc(Int,(v_["A"]/v_["B"]))
            v_["B/A"] = trunc(Int,(v_["B"]/v_["A"]))
            v_["SUM"] = v_["A/B"]+v_["B/A"]
            v_["A"] = v_["A"]*v_["A/B"]
            v_["B"] = v_["B"]*v_["B/A"]
            v_["MAXA,B"] = v_["A"]+v_["B"]
            v_["MAXA,B"] = trunc(Int,(v_["MAXA,B"]/v_["SUM"]))
            v_["MAXA,B"] = v_["MAXA,B"]-v_["ABSA+B"]
            v_["MAXI-5,1"] = v_["MAXA,B"]
            v_["A"] = v_["I+1"]
            v_["B"] = v_["N"]
            v_["A"] = -1*v_["A"]
            v_["B"] = -1*v_["B"]
            v_["A"] = Float64(v_["A"])
            v_["ABSA"] = abs(v_["A"])
            v_["ABSA"] = trunc(Int,v_["ABSA"])
            v_["B"] = Float64(v_["B"])
            v_["ABSB"] = abs(v_["B"])
            v_["ABSB"] = trunc(Int,v_["ABSB"])
            v_["ABSA+B"] = v_["ABSA"]+v_["ABSB"]
            v_["A"] = v_["A"]+v_["ABSA+B"]
            v_["B"] = v_["B"]+v_["ABSA+B"]
            v_["A/B"] = trunc(Int,(v_["A"]/v_["B"]))
            v_["B/A"] = trunc(Int,(v_["B"]/v_["A"]))
            v_["SUM"] = v_["A/B"]+v_["B/A"]
            v_["A"] = v_["A"]*v_["A/B"]
            v_["B"] = v_["B"]*v_["B/A"]
            v_["MAXA,B"] = v_["A"]+v_["B"]
            v_["MAXA,B"] = trunc(Int,(v_["MAXA,B"]/v_["SUM"]))
            v_["MINA,B"] = v_["ABSA+B"]-v_["MAXA,B"]
            v_["MINI+1,N"] = v_["MINA,B"]
            for J = Int64(v_["MAXI-5,1"]):Int64(v_["MINI+1,N"])
                ig = ig_["OBJ"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["S"*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(1.0))
            end
        end
        for K = Int64(v_["1"]):Int64(v_["N/2"])
            ig = ig_["C"*string(K)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["P"*string(K)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               6.26382E+04
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AY-V-V"
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

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]
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

    elseif action == "eCUBE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]^3
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*EV_[1]^2
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

    elseif action == "eXEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        U_[2,1] = U_[2,1]+1
        U_[2,2] = U_[2,2]-1
        U_[2,3] = U_[2,3]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        EXPW = exp(IV_[2])
        UEXPW = IV_[1]*EXPW
        f_   = UEXPW
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPW
            g_[2] = UEXPW
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = EXPW
                H_[2,1] = H_[1,2]
                H_[2,2] = UEXPW
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

    elseif action == "gL7d3"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        Z = abs(GVAR_)
        f_= Z^(7.0/3.0)
        if nargout>1
            g_ = 7.0*sign(GVAR_)*(Z^(4.0/3.0))/3.0
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 28.0*Z^(1.0/3.0)/9.0
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

