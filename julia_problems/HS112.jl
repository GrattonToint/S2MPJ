function HS112(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS112
#    *********
# 
#    This problem is a chemical equilibrium problem involving 3 linear
#    equality constraints.
# 
#    Source: problem 80 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: A.R. Conn, Mar 1990.
# 
#    classification = "C-OLR2-MY-10-3"
# 
#    N is the number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS112"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 10
        v_["1"] = 1
        v_["C1"] = -6.089
        v_["C2"] = -17.164
        v_["C3"] = -34.054
        v_["C4"] = -5.914
        v_["C5"] = -24.721
        v_["C6"] = -14.986
        v_["C7"] = -24.100
        v_["C8"] = -10.708
        v_["C9"] = -26.662
        v_["C10"] = -22.179
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
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["C"*string(I)])
        end
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("CON2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON2")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("CON3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON3")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
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
        pbm.gconst[ig_["CON1"]] = Float64(2.0)
        pbm.gconst[ig_["CON2"]] = Float64(1.0)
        pbm.gconst[ig_["CON3"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.0e-6,pb.n)
        pb.xupper = fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.1),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eLOG", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eLOGSUM", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"X1")
        loaset(elftv,it,3,"X2")
        loaset(elftv,it,4,"X3")
        loaset(elftv,it,5,"X4")
        loaset(elftv,it,6,"X5")
        loaset(elftv,it,7,"X6")
        loaset(elftv,it,8,"X7")
        loaset(elftv,it,9,"X8")
        loaset(elftv,it,10,"X9")
        loaset(elftv,it,11,"X10")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "XLOGX"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOG")
            arrset(ielftype,ie,iet_["eLOG"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "XLOGS"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOGSUM")
            arrset(ielftype,ie,iet_["eLOGSUM"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X7",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X8",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X9"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X9",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X10"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-6),nothing,Float64(0.1))
            posev = findfirst(x->x=="X10",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XLOGX"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["XLOGS"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -47.707579
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
        pb.pbclass = "C-OLR2-MY-10-3"
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

    elseif action == "eLOG"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LOGX = log(EV_[1])
        f_   = EV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = LOGX+1.0
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.0/EV_[1]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eLOGSUM"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,11)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[2,2] = U_[2,2]+1
        U_[2,3] = U_[2,3]+1
        U_[2,4] = U_[2,4]+1
        U_[2,5] = U_[2,5]+1
        U_[2,6] = U_[2,6]+1
        U_[2,7] = U_[2,7]+1
        U_[2,8] = U_[2,8]+1
        U_[2,9] = U_[2,9]+1
        U_[2,10] = U_[2,10]+1
        U_[2,11] = U_[2,11]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        LOGSUM = log(IV_[2])
        f_   = -IV_[1]*LOGSUM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -LOGSUM
            g_[2] = -IV_[1]/IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -1.0/IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = IV_[1]/IV_[2]^2
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

