function LEWISPOL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Adrian Lewis Polynomial Problem,
#    The problem is a transformation of a number theory integer
#    programming problem.
# 
#    Source:
#    A. Lewis, private communication.
# 
#    SIF input: A.R. Conn and Ph. Toint, March 1990.
# 
#    classification = "C-QOR2-AN-6-9"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LEWISPOL"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 6
        v_["DEG"] = 3
        v_["PEN"] = 1.0e4
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["DEG-1"] = -1+v_["DEG"]
        v_["N-1"] = -1+v_["N"]
        v_["N+1"] = 1+v_["N"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            iv,ix_,_ = s2mpj_ii("A"*string(J),ix_)
            arrset(pb.xnames,iv,"A"*string(J))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            v_["C"*string(Int64(v_["0"]))*","*string(J)] = 1.0
            ig,ig_,_ = s2mpj_ii("D0",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"D0")
            iv = ix_["A"*string(J)]
            pbm.A[ig,iv] += Float64(v_["C"*string(Int64(v_["0"]))*","*string(J)])
        end
        for I = Int64(v_["1"]):Int64(v_["DEG-1"])
            v_["I-1"] = -1+I
            for J = Int64(I):Int64(v_["N-1"])
                v_["RJ"] = Float64(J)
                v_["C"*string(I)*","*string(J)]  = (
                      v_["C"*string(Int64(v_["I-1"]))*","*string(J)]*v_["RJ"])
                ig,ig_,_ = s2mpj_ii("D"*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"D"*string(I))
                iv = ix_["A"*string(J)]
                pbm.A[ig,iv] += Float64(v_["C"*string(I)*","*string(J)])
            end
        end
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("INT"*string(J),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"INT"*string(J))
            iv = ix_["A"*string(J)]
            pbm.A[ig,iv] += Float64(-1.0)
            arrset(pbm.gscale,ig,Float64(v_["PEN"]))
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
        v_["CT"*string(Int64(v_["0"]))] = -1.0
        pbm.gconst[ig_["D0"]] = Float64(v_["CT"*string(Int64(v_["0"]))])
        for I = Int64(v_["1"]):Int64(v_["DEG-1"])
            v_["I-1"] = -1+I
            v_["-I"] = -1*I
            v_["N+1-I"] = v_["N+1"]+v_["-I"]
            v_["VAL"] = Float64(v_["N+1-I"])
            v_["CT"*string(I)] = v_["CT"*string(Int64(v_["I-1"]))]*v_["VAL"]
            pbm.gconst[ig_["D"*string(I)]] = Float64(v_["CT"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(-10.0,pb.n)
        pb.xupper = fill(10.0,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"A0")
            pb.x0[ix_["A0"]] = Float64(-1.0)
        else
            pb.y0[findfirst(x->x==ig_["A0"],pbm.congrps)] = Float64(-1.0)
        end
        if haskey(ix_,"A1")
            pb.x0[ix_["A1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["A1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"A2")
            pb.x0[ix_["A2"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["A2"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"A3")
            pb.x0[ix_["A3"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["A3"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"A4")
            pb.x0[ix_["A4"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["A4"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"A5")
            pb.x0[ix_["A5"]] = Float64(-1.0)
        else
            pb.y0[findfirst(x->x==ig_["A5"],pbm.congrps)] = Float64(-1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eCB", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            ename = "O"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "A"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-10.0),Float64(10.0),nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCB")
            arrset(ielftype,ie,iet_["eCB"])
            vname = "A"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(-10.0),Float64(10.0),nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for J = Int64(v_["0"]):Int64(v_["N-1"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["O"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            ig = ig_["INT"*string(J)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
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
        pb.pbclass = "C-QOR2-AN-6-9"
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

    elseif action == "eCB"

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

