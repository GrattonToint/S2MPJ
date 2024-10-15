function HIMMELBK(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELBK
#    *********
# 
#    A blending problem for multi-component mixtures, by Paviani.
#    It has a linear objective and linear and nonlinear constraints.
# 
#    Compared to the problem specified in Himmelblau, the inequality
#    constraints have been removed, because, as stated in this source,
#    they impose that
#    X(1)=X(2)=X(3)=X(7)=X(9)=X(9)=X(13)=X(14)=X(15)=X(19)=X(20)=X(21)=0
#    which is clearly contradictory with the given solution(s).  As
#    there does not seem to be a natural way to correct this statement
#    without knowing more about the original problem, the troublesome
#    constraints have been removed.
# 
#    Source: from problem 20 in
#    D.H. Himmelblau,
#    "Applied nonlinear programming",
#    McGraw-Hill, New-York, 1972.
# 
#    SIF input: Ph. Toint, March 1991.
# 
#    classification = "C-LOR2-MN-24-14"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HIMMELBK"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["F"] = 142.22471
        v_["B1"] = 44.094
        v_["B2"] = 58.12
        v_["B3"] = 58.12
        v_["B4"] = 137.4
        v_["B5"] = 120.9
        v_["B6"] = 170.9
        v_["B7"] = 62.501
        v_["B8"] = 84.94
        v_["B9"] = 133.425
        v_["B10"] = 82.507
        v_["B11"] = 46.07
        v_["B12"] = 60.097
        v_["B13"] = 44.094
        v_["B14"] = 58.12
        v_["B15"] = 58.12
        v_["B16"] = 137.4
        v_["B17"] = 120.9
        v_["B18"] = 170.9
        v_["B19"] = 62.501
        v_["B20"] = 84.94
        v_["B21"] = 133.425
        v_["B22"] = 82.507
        v_["B23"] = 46.07
        v_["B24"] = 60.097
        v_["C1"] = 123.7
        v_["C2"] = 31.7
        v_["C3"] = 45.7
        v_["C4"] = 14.7
        v_["C5"] = 84.7
        v_["C6"] = 27.7
        v_["C7"] = 49.7
        v_["C8"] = 7.1
        v_["C9"] = 2.1
        v_["C10"] = 17.7
        v_["C11"] = 0.85
        v_["C12"] = 0.64
        v_["D1"] = 123.7
        v_["D2"] = 31.7
        v_["D3"] = 45.7
        v_["D4"] = 14.7
        v_["D5"] = 84.7
        v_["D6"] = 27.7
        v_["D7"] = 49.7
        v_["D8"] = 7.1
        v_["D9"] = 2.1
        v_["D10"] = 17.7
        v_["D11"] = 0.85
        v_["D12"] = 0.64
        v_["1"] = 1
        v_["12"] = 12
        v_["13"] = 13
        v_["24"] = 24
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for K = Int64(v_["1"]):Int64(v_["24"])
            iv,ix_,_ = s2mpj_ii("X"*string(K),ix_)
            arrset(pb.xnames,iv,"X"*string(K))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(0.0693)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(0.0577)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(0.05)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.2)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(0.26)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(0.55)
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(0.06)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(0.1)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(0.12)
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(0.18)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(0.1)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(0.09)
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(0.0693)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(0.0577)
        iv = ix_["X15"]
        pbm.A[ig,iv] += Float64(0.05)
        iv = ix_["X16"]
        pbm.A[ig,iv] += Float64(0.2)
        iv = ix_["X17"]
        pbm.A[ig,iv] += Float64(0.26)
        iv = ix_["X18"]
        pbm.A[ig,iv] += Float64(0.55)
        iv = ix_["X19"]
        pbm.A[ig,iv] += Float64(0.06)
        iv = ix_["X20"]
        pbm.A[ig,iv] += Float64(0.1)
        iv = ix_["X21"]
        pbm.A[ig,iv] += Float64(0.12)
        iv = ix_["X22"]
        pbm.A[ig,iv] += Float64(0.18)
        iv = ix_["X23"]
        pbm.A[ig,iv] += Float64(0.1)
        iv = ix_["X24"]
        pbm.A[ig,iv] += Float64(0.09)
        for I = Int64(v_["1"]):Int64(v_["12"])
            ig,ig_,_ = s2mpj_ii("CA"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CA"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["24"])
            ig,ig_,_ = s2mpj_ii("CA13",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CA13")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I+12"] = 12+I
            v_["1/DI"] = 1.0/v_["D"*string(I)]
            ig,ig_,_ = s2mpj_ii("CA14",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"CA14")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/DI"])
            v_["F/BI+12"] = v_["F"]/v_["B"*string(Int64(v_["I+12"]))]
            iv = ix_["X"*string(Int64(v_["I+12"]))]
            pbm.A[ig,iv] += Float64(v_["F/BI+12"])
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
        pbm.gconst[ig_["CA13"]] = Float64(1.0)
        pbm.gconst[ig_["CA14"]] = Float64(1.671)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.04),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I+12"] = 12+I
            for J = Int64(v_["1"]):Int64(v_["12"])
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
                vname = "X"*string(Int64(v_["I+12"]))
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.04))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.04))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
            for J = Int64(v_["13"]):Int64(v_["24"])
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"en2PR")
                arrset(ielftype,ie,iet_["en2PR"])
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.04))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.04))
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["12"])
            v_["I+12"] = 12+I
            for J = Int64(v_["1"]):Int64(v_["12"])
                v_["BI/BJ"] = v_["B"*string(I)]/v_["B"*string(J)]
                v_["40BI/BJ"] = 40.0*v_["BI/BJ"]
                ig = ig_["CA"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["40BI/BJ"]))
            end
            for J = Int64(v_["13"]):Int64(v_["24"])
                v_["B+/BJ"] = v_["B"*string(Int64(v_["I+12"]))]/v_["B"*string(J)]
                v_["CB+/BJ"] = v_["C"*string(I)]*v_["B+/BJ"]
                v_["-CB+/BJ"] = -1.0*v_["CB+/BJ"]
                ig = ig_["CA"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(v_["-CB+/BJ"]))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                0.0893344
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-LOR2-MN-24-14"
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

