function HS118(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS118
#    *********
# 
#    Source: problem 118 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: B Baudson, Jan 1990.
# 
#    classification = "C-QLR2-AN-15-17"
# 
#    Other useful parameters
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS118"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["0"] = 0
        v_["1"] = 1
        v_["4"] = 4
        v_["8"] = 8
        v_["12"] = 12
        v_["15"] = 15
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["15"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for K = Int64(v_["0"]):Int64(v_["4"])
            v_["3K"] = 3*K
            v_["3K+1"] = 1+v_["3K"]
            v_["3K+2"] = 2+v_["3K"]
            v_["3K+3"] = 3+v_["3K"]
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["3K+1"]))]
            pbm.A[ig,iv] += Float64(2.3)
            iv = ix_["X"*string(Int64(v_["3K+2"]))]
            pbm.A[ig,iv] += Float64(1.7)
            iv = ix_["X"*string(Int64(v_["3K+3"]))]
            pbm.A[ig,iv] += Float64(2.2)
        end
        for K = Int64(v_["1"]):Int64(v_["4"])
            v_["3K"] = 3*K
            v_["3K+1"] = 1+v_["3K"]
            v_["3K+2"] = 2+v_["3K"]
            v_["3K+3"] = 3+v_["3K"]
            v_["3K-2"] = -2+v_["3K"]
            v_["3K-1"] = -1+v_["3K"]
            ig,ig_,_ = s2mpj_ii("A"*string(K),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"A"*string(K))
            iv = ix_["X"*string(Int64(v_["3K+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["3K-2"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("B"*string(K),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"B"*string(K))
            iv = ix_["X"*string(Int64(v_["3K+3"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["3K"]))]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(K),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"C"*string(K))
            iv = ix_["X"*string(Int64(v_["3K+2"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(Int64(v_["3K-1"]))]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        ig,ig_,_ = s2mpj_ii("D1",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"D1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("D2",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"D2")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("D3",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"D3")
        iv = ix_["X7"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X8"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X9"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("D4",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"D4")
        iv = ix_["X10"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X11"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X12"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("D5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"D5")
        iv = ix_["X13"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X14"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X15"]
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
        for K = Int64(v_["1"]):Int64(v_["4"])
            pbm.gconst[ig_["A"*string(K)]] = Float64(-7.0)
            pbm.gconst[ig_["B"*string(K)]] = Float64(-7.0)
            pbm.gconst[ig_["C"*string(K)]] = Float64(-7.0)
        end
        pbm.gconst[ig_["D1"]] = Float64(60.0)
        pbm.gconst[ig_["D2"]] = Float64(50.0)
        pbm.gconst[ig_["D3"]] = Float64(70.0)
        pbm.gconst[ig_["D4"]] = Float64(85.0)
        pbm.gconst[ig_["D5"]] = Float64(100.0)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for K = Int64(v_["1"]):Int64(v_["4"])
            arrset(grange,ig_["A"*string(K)],Float64(13.0))
            arrset(grange,ig_["B"*string(K)],Float64(13.0))
            arrset(grange,ig_["C"*string(K)],Float64(14.0))
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.0,pb.n)
        pb.xupper = fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 8.0
        pb.xupper[ix_["X1"]] = 21.0
        pb.xlower[ix_["X2"]] = 43.0
        pb.xupper[ix_["X2"]] = 57.0
        pb.xlower[ix_["X3"]] = 3.0
        pb.xupper[ix_["X3"]] = 16.0
        for K = Int64(v_["1"]):Int64(v_["4"])
            v_["3K"] = 3*K
            v_["3K+1"] = 1+v_["3K"]
            v_["3K+2"] = 2+v_["3K"]
            v_["3K+3"] = 3+v_["3K"]
            pb.xupper[ix_["X"*string(Int64(v_["3K+1"]))]] = 90.0
            pb.xupper[ix_["X"*string(Int64(v_["3K+2"]))]] = 120.0
            pb.xupper[ix_["X"*string(Int64(v_["3K+3"]))]] = 60.0
        end
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(20.0),pb.n)
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(55.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(55.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(15.0)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(15.0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(60.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(60.0)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(60.0)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(60.0)
        end
        if haskey(ix_,"X11")
            pb.x0[ix_["X11"]] = Float64(60.0)
        else
            pb.y0[findfirst(x->x==ig_["X11"],pbm.congrps)] = Float64(60.0)
        end
        if haskey(ix_,"X14")
            pb.x0[ix_["X14"]] = Float64(60.0)
        else
            pb.y0[findfirst(x->x==ig_["X14"],pbm.congrps)] = Float64(60.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["15"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),nothing,Float64(20.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for K = Int64(v_["0"]):Int64(v_["4"])
            v_["3K"] = 3*K
            v_["3K+1"] = 1+v_["3K"]
            v_["3K+2"] = 2+v_["3K"]
            v_["3K+3"] = 3+v_["3K"]
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3K+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(0.0001))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3K+2"]))])
            loaset(pbm.grelw,ig,posel,Float64(0.0001))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["3K+3"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(0.00015))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               664.82045
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QLR2-AN-15-17"
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

    elseif action == "eSQ"

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

