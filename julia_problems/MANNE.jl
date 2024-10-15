function MANNE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MANNE
#    *********
# 
#    A variable dimension econometric equilibrium problem
#    suggested by A. Manne
# 
#    Source:
#    B. Murtagh and M. Saunders,
#    Mathematical Programming Studies 16, pp. 84-117,
#    (example 5.12).
# 
#    SIF input: N. Gould and Ph. Toint, March 1990.
# 
#    classification = "C-OOR2-MN-V-V"
# 
#    Number of periods
#    The number of variables in the problem N = 3*T
# 
#       Alternative values for the SIF file parameters:
# IE T                   100            $-PARAMETER n = 300    original value
# IE T                   365            $-PARAMETER n = 995
# IE T                   1000           $-PARAMETER n = 3000
# IE T                   2000           $-PARAMETER n = 6000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MANNE"

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
            v_["T"] = Int64(4);  #  SIF file default value
        else
            v_["T"] = Int64(args[1]);
        end
        v_["GROW"] = 0.03
        v_["BETA"] = 0.95
        v_["XK0"] = 3.0
        v_["XC0"] = 0.95
        v_["XI0"] = 0.05
        v_["B"] = 0.25
        v_["BPROB"] = 1.0
        v_["1"] = 1
        v_["2"] = 2
        v_["T-1"] = -1+v_["T"]
        v_["T-2"] = -2+v_["T"]
        v_["LOGXK"] = log(v_["XK0"])
        v_["BLOGX"] = v_["LOGXK"]*v_["B"]
        v_["XK0**B"] = exp(v_["BLOGX"])
        v_["NUM"] = v_["XC0"]+v_["XI0"]
        v_["A"] = v_["NUM"]/v_["XK0**B"]
        v_["1-B"] = 1.0-v_["B"]
        v_["1+G"] = 1.0+v_["GROW"]
        v_["LOG1+G"] = log(v_["1+G"])
        v_["SOME"] = v_["LOG1+G"]*v_["1-B"]
        v_["GFAC"] = exp(v_["SOME"])
        v_["AT1"] = v_["A"]*v_["GFAC"]
        v_["BT1"] = 0.0+v_["BETA"]
        for J = Int64(v_["2"]):Int64(v_["T"])
            v_["J-1"] = -1+J
            v_["AT"*string(J)] = v_["AT"*string(Int64(v_["J-1"]))]*v_["GFAC"]
            v_["BT"*string(J)] = v_["BT"*string(Int64(v_["J-1"]))]*v_["BETA"]
        end
        v_["1-BETA"] = 1.0-v_["BETA"]
        v_["1/1-BETA"] = 1.0/v_["1-BETA"]
        v_["BT"*string(Int64(v_["T"]))]  = (
              v_["BT"*string(Int64(v_["T"]))]*v_["1/1-BETA"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["T"])
            iv,ix_,_ = s2mpj_ii("C"*string(I),ix_)
            arrset(pb.xnames,iv,"C"*string(I))
            iv,ix_,_ = s2mpj_ii("I"*string(I),ix_)
            arrset(pb.xnames,iv,"I"*string(I))
            iv,ix_,_ = s2mpj_ii("K"*string(I),ix_)
            arrset(pb.xnames,iv,"K"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["T"])
            ig,ig_,_ = s2mpj_ii("NL"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"NL"*string(I))
            iv = ix_["C"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["I"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["T-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("L"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"L"*string(I))
            iv = ix_["K"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["K"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["I"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        ig,ig_,_ = s2mpj_ii("L"*string(Int64(v_["T"])),ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L"*string(Int64(v_["T"])))
        iv = ix_["K"*string(Int64(v_["T"]))]
        pbm.A[ig,iv] += Float64(v_["GROW"])
        ig,ig_,_ = s2mpj_ii("L"*string(Int64(v_["T"])),ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L"*string(Int64(v_["T"])))
        iv = ix_["I"*string(Int64(v_["T"]))]
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
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        grange[gegrps,1] = fill(Inf,pb.nge)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["K1"]] = 3.05
        pb.xupper[ix_["K1"]] = 3.05
        for I = Int64(v_["2"]):Int64(v_["T"])
            pb.xlower[ix_["K"*string(I)]] = 3.05
        end
        v_["1.04**T"] = 0.05
        for I = Int64(v_["1"]):Int64(v_["T"])
            v_["1.04**T"] = 1.04*v_["1.04**T"]
            pb.xlower[ix_["C"*string(I)]] = 0.95
            pb.xlower[ix_["I"*string(I)]] = 0.05
            pb.xupper[ix_["I"*string(I)]] = v_["1.04**T"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"K1")
            pb.x0[ix_["K1"]] = Float64(3.05)
        else
            pb.y0[findfirst(x->x==ig_["K1"],pbm.congrps)] = Float64(3.05)
        end
        for I = Int64(v_["2"]):Int64(v_["T"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["I-1/10"] = 0.1*v_["RI-1"]
            v_["VAL"] = 3.0+v_["I-1/10"]
            if haskey(ix_,"K"*string(I))
                pb.x0[ix_["K"*string(I)]] = Float64(v_["VAL"])
            else
                pb.y0[findfirst(x->x==ig_["K"*string(I)],pbm.congrps)] = Float64(v_["VAL"])
            end
        end
        for I = Int64(v_["1"]):Int64(v_["T"])
            if haskey(ix_,"C"*string(I))
                pb.x0[ix_["C"*string(I)]] = Float64(0.95)
            else
                pb.y0[findfirst(x->x==ig_["C"*string(I)],pbm.congrps)] = Float64(0.95)
            end
            if haskey(ix_,"I"*string(I))
                pb.x0[ix_["I"*string(I)]] = Float64(0.05)
            else
                pb.y0[findfirst(x->x==ig_["I"*string(I)],pbm.congrps)] = Float64(0.05)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eLOGS", iet_)
        loaset(elftv,it,1,"C")
        it,iet_,_ = s2mpj_ii( "ePOWER", iet_)
        loaset(elftv,it,1,"K")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"B")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["T"])
            ename = "LOGC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOGS")
            arrset(ielftype,ie,iet_["eLOGS"])
            vname = "C"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="C",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "KS"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePOWER")
            arrset(ielftype,ie,iet_["ePOWER"])
            vname = "K"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="K",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="B",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["T"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["LOGC"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["BT"*string(I)]))
            ig = ig_["NL"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["KS"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["AT"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -9.7457259D-01
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-MN-V-V"
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

    elseif action == "eLOGS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = log(EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -1.0/EV_[1]^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePOWER"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]^pbm.elpar[iel_][1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EV_[1]^(pbm.elpar[iel_][1]-1.0)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1]  = (
                      pbm.elpar[iel_][1]*(pbm.elpar[iel_][1]-1.0)*EV_[1]^(pbm.elpar[iel_][1]-2.0))
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

