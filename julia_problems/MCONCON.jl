function MCONCON(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Another small gas network problem.
# 
#    SIF input: Sybille Schachler, Oxford, August 1992.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-LOI2-MN-15-11"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MCONCON"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 7
        v_["M"] = 4
        v_["M+1"] = 1+v_["M"]
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
            iv,ix_,_ = s2mpj_ii("Q"*string(I),ix_)
            arrset(pb.xnames,iv,"Q"*string(I))
            iv,ix_,_ = s2mpj_ii("F"*string(I),ix_)
            arrset(pb.xnames,iv,"F"*string(I))
        end
        for I = Int64(v_["M+1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("OBJECT",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["P"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("PAN"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PAN"*string(I))
        end
        ig,ig_,_ = s2mpj_ii("MBAL1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL1")
        iv = ix_["Q1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["F3"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("MBAL2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL2")
        iv = ix_["Q1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["F1"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("MBAL3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL3")
        iv = ix_["Q2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["F1"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("MBAL4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL4")
        iv = ix_["Q2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["Q3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("MBAL5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL5")
        iv = ix_["Q3"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["F2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("MBAL6",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL6")
        iv = ix_["Q4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["F2"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("MBAL7",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"MBAL7")
        iv = ix_["Q4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["F4"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        v_["DEMAND"] = -1000.0
        pbm.gconst[ig_["MBAL4"]] = Float64(v_["DEMAND"])
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        v_["PMAX1"] = 914.73
        v_["PMAX2"] = 904.73
        pb.xupper[ix_["P3"]] = v_["PMAX2"]
        pb.xupper[ix_["P5"]] = v_["PMAX2"]
        pb.xupper[ix_["P1"]] = v_["PMAX1"]
        pb.xupper[ix_["P7"]] = v_["PMAX1"]
        pb.xupper[ix_["F4"]] = 400.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["P"*string(I)]] = Float64(965.0)
        end
        pb.x0[ix_["Q1"]] = Float64(100.0)
        pb.x0[ix_["Q2"]] = Float64(100.0)
        pb.x0[ix_["Q3"]] = Float64(-100.0)
        pb.x0[ix_["Q4"]] = Float64(-100.0)
        pb.x0[ix_["F1"]] = Float64(1000.0)
        pb.x0[ix_["F2"]] = Float64(1000.0)
        pb.x0[ix_["F3"]] = Float64(1000.0)
        pb.x0[ix_["F4"]] = Float64(1000.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eFORQ", iet_)
        loaset(elftv,it,1,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "PSQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQR")
            arrset(ielftype,ie,iet_["eSQR"])
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "QTO"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eFORQ")
            arrset(ielftype,ie,iet_["eFORQ"])
            vname = "Q"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        v_["K"] = -0.597053452
        ig = ig_["PAN1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ2"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["QTO1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["K"]))
        ig = ig_["PAN2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ4"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["QTO2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["K"]))
        ig = ig_["PAN3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ5"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["QTO3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["K"]))
        ig = ig_["PAN4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ6"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["PSQ7"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["QTO4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["K"]))
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
        pb.pbclass = "C-LOI2-MN-15-11"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*abs(EV_[1])
        PO = EV_[1]>0.0
        if PO
            GO = 2*EV_[1]
        end
        if !PO
            GO = -2*EV_[1]
        end
        if PO
            HO = 2
        end
        if !PO
            HO = -2
        end
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = GO
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = HO
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eFORQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*abs(EV_[1])^0.8539
        POS = EV_[1]>0.0
        if POS
            GG = 1.8539*EV_[1]^0.8539
        end
        if !POS
            GG = 1.8539*abs(EV_[1])^0.8539
        end
        if POS
            HH = 1.8539*0.8539*EV_[1]^(-0.1461)
        end
        if !POS
            HH = -1.8539*0.8539*abs(EV_[1])^(-0.1461)
        end
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = GG
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = HH
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

