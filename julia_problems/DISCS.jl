function DISCS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DISCS
#    *********
# 
#    The problem is to replace the nodes of a planar graph by discs
#    such that adjacent nodes have touching discs and non adjacent nodes
#    disjoint discs.  The smallest disc is constrained to have radius equal to
#    one and is centered at the origin.  The second discs is also
#    constrained to have its Y = 0.0, in order to cancel invariance under
#    rotation.
# 
#    Source:
#    W. Pulleyblank,
#    private communication, 1991.
# 
#    classification = "C-LQR2-MY-36-66"
# 
#    SIF input: A.R. Conn and Ph. Toint, April 1991.
# 
#    Number of nodes
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DISCS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["NNODES"] = 12
        v_["EPSIL"] = 0.0001
        v_["1"] = 1
        v_["2"] = 2
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["A"*string(J)*","*string(I)] = 0.0
            end
        end
        v_["A1,2"] = 1.01
        v_["A1,7"] = 1.01
        v_["A2,3"] = 1.01
        v_["A2,4"] = 1.01
        v_["A3,4"] = 1.01
        v_["A4,5"] = 1.01
        v_["A4,6"] = 1.01
        v_["A5,6"] = 1.01
        v_["A5,11"] = 1.01
        v_["A6,7"] = 1.01
        v_["A7,8"] = 1.01
        v_["A7,9"] = 1.01
        v_["A8,9"] = 1.01
        v_["A8,10"] = 1.01
        v_["A9,10"] = 1.01
        v_["A10,11"] = 1.01
        v_["A10,12"] = 1.01
        v_["A11,12"] = 1.01
        v_["RNODES"] = Float64(v_["NNODES"])
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2mpj_ii("R"*string(I),ix_)
            arrset(pb.xnames,iv,"R"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            ig,ig_,_ = s2mpj_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["2"]):Int64(v_["NNODES"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["RAIJ"] = v_["A"*string(J)*","*string(I)]
                v_["AIJ"] = trunc(Int,v_["RAIJ"])
                for K = Int64(v_["1"]):Int64(v_["AIJ"])
                    ig,ig_,_ = s2mpj_ii("S"*string(J)*","*string(I),ig_)
                    arrset(gtype,ig,"==")
                    arrset(pb.cnames,ig,"S"*string(J)*","*string(I))
                end
                v_["AIJ-1"] = -1+v_["AIJ"]
                v_["NAIJ"] = -1*v_["AIJ-1"]
                for K = Int64(v_["1"]):Int64(v_["NAIJ"])
                    ig,ig_,_ = s2mpj_ii("S"*string(J)*","*string(I),ig_)
                    arrset(gtype,ig,"<=")
                    arrset(pb.cnames,ig,"S"*string(J)*","*string(I))
                end
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
        v_["-EPSIL"] = -1.0*v_["EPSIL"]
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                v_["RAIJ"] = v_["A"*string(J)*","*string(I)]
                v_["AIJ"] = trunc(Int,v_["RAIJ"])
                v_["AIJ-1"] = -1+v_["AIJ"]
                v_["NAIJ"] = -1*v_["AIJ-1"]
                for K = Int64(v_["1"]):Int64(v_["NAIJ"])
                    pbm.gconst[ig_["S"*string(J)*","*string(I)]] = Float64(v_["-EPSIL"])
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            pb.xlower[ix_["R"*string(I)]] = 1.0
        end
        pb.xlower[ix_["X"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["1"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["1"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["2"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["2"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["NNODES"])
            v_["RI"] = Float64(I)
            v_["RM"] = 0.03*v_["RI"]
            v_["RP"] = 0.0055555*v_["RI"]
            pb.x0[ix_["R"*string(I)]] = Float64(v_["RI"])
            pb.x0[ix_["X"*string(I)]] = Float64(v_["RM"])
            pb.x0[ix_["Y"*string(I)]] = Float64(v_["RP"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eIPSQ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eIMSQ", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["NNODES"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ename = "DR"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eIPSQ")
                arrset(ielftype,ie,iet_["eIPSQ"])
                vname = "R"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "R"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DX"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eIMSQ")
                arrset(ielftype,ie,iet_["eIMSQ"])
                vname = "X"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "DY"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eIMSQ")
                arrset(ielftype,ie,iet_["eIMSQ"])
                vname = "Y"*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["NNODES"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["I-1"])
                ig = ig_["S"*string(J)*","*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["DR"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["DX"*string(J)*","*string(I)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["DY"*string(J)*","*string(I)])
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# ZL DISCS                              RNODES
#    Solution
# LO SOLTN(12)           20.46122911
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LQR2-MY-36-66"
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

    elseif action == "eIPSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]+1
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

    elseif action == "eIMSQ"

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

