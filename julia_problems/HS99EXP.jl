function HS99EXP(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS99EXP
#    *********
# 
#    Source: an expanded form of problem 99 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-OOR2-AN-31-21"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS99EXP"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["T1"] = 0.0
        v_["T2"] = 25.0
        v_["T3"] = 50.0
        v_["T4"] = 100.0
        v_["T5"] = 150.0
        v_["T6"] = 200.0
        v_["T7"] = 290.0
        v_["T8"] = 380.0
        v_["A1"] = 0.0
        v_["A2"] = 50.0
        v_["A3"] = 50.0
        v_["A4"] = 75.0
        v_["A5"] = 75.0
        v_["A6"] = 75.0
        v_["A7"] = 100.0
        v_["A8"] = 100.0
        v_["B"] = 32.0
        v_["1"] = 1
        v_["2"] = 2
        v_["7"] = 7
        v_["8"] = 8
        for I = Int64(v_["2"]):Int64(v_["8"])
            v_["I-1"] = -1+I
            v_["DT"*string(I)] = v_["T"*string(I)]-v_["T"*string(Int64(v_["I-1"]))]
            v_["DTISQ"] = v_["DT"*string(I)]*v_["DT"*string(I)]
            v_["DT"*string(I)] = 0.5*v_["DTISQ"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["7"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("R"*string(I),ix_)
            arrset(pb.xnames,iv,"R"*string(I))
            iv,ix_,_ = s2mpj_ii("Q"*string(I),ix_)
            arrset(pb.xnames,iv,"Q"*string(I))
            iv,ix_,_ = s2mpj_ii("S"*string(I),ix_)
            arrset(pb.xnames,iv,"S"*string(I))
        end
        iv,ix_,_ = s2mpj_ii("R"*string(Int64(v_["8"])),ix_)
        arrset(pb.xnames,iv,"R"*string(Int64(v_["8"])))
        iv,ix_,_ = s2mpj_ii("Q"*string(Int64(v_["8"])),ix_)
        arrset(pb.xnames,iv,"Q"*string(Int64(v_["8"])))
        iv,ix_,_ = s2mpj_ii("S"*string(Int64(v_["8"])),ix_)
        arrset(pb.xnames,iv,"S"*string(Int64(v_["8"])))
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["2"]):Int64(v_["8"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("R"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"R"*string(I))
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["R"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("Q"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"Q"*string(I))
            iv = ix_["Q"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["Q"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["S"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["DT"*string(I)])
            ig,ig_,_ = s2mpj_ii("S"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"S"*string(I))
            iv = ix_["S"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["S"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(1.0)
        end
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["R"*string(Int64(v_["8"]))]
        pbm.A[ig,iv] += Float64(1.0)
        arrset(pbm.gscale,ig,Float64(-1.0))
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
        for I = Int64(v_["2"]):Int64(v_["7"])
            v_["RHS"] = v_["DT"*string(I)]*v_["B"]
            pbm.gconst[ig_["Q"*string(I)]] = Float64(v_["RHS"])
            v_["RHS"] = v_["DT"*string(I)]*v_["B"]
            pbm.gconst[ig_["S"*string(I)]] = Float64(v_["RHS"])
        end
        pbm.gconst[ig_["Q"*string(Int64(v_["8"]))]] = Float64(100000.0)
        pbm.gconst[ig_["S"*string(Int64(v_["8"]))]] = Float64(1000.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["R"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["R"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Q"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Q"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["S"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["S"*string(Int64(v_["1"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["7"])
            pb.xlower[ix_["X"*string(I)]] = 0.0
            pb.xupper[ix_["X"*string(I)]] = 1.58
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["7"])
            pb.x0[ix_["X"*string(I)]] = Float64(0.5)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSN", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eCS", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["7"])
            ename = "SNX"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSN")
            arrset(ielftype,ie,iet_["eSN"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "CSX"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCS")
            arrset(ielftype,ie,iet_["eCS"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        arrset(pbm.grftype,ig,"gL2")
        for I = Int64(v_["2"]):Int64(v_["8"])
            v_["I-1"] = -1+I
            v_["W"] = v_["A"*string(I)]*v_["DT"*string(I)]
            ig = ig_["R"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CSX"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["W"]))
            ig = ig_["S"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["SNX"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["W"]))
            v_["W"] = v_["A"*string(I)]*v_["DT"*string(I)]
            ig = ig_["Q"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["SNX"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["W"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -831079892.0
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
        pb.pbclass = "C-OOR2-AN-31-21"
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

    elseif action == "eSN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SNX = sin(EV_[1])
        f_   = SNX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -SNX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CSX = cos(EV_[1])
        f_   = CSX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -sin(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -CSX
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

