function OET6(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OET6
#    *********
# 
#    A nonlinear programming formulation of a discretization of
#    a nonlinear Chebychev problem.
# 
#    The problem is
# 
#        min  max || phi(x,w) ||, for all w in the interval I.
#         x    w
# 
#    I is discretized, and the problem solved over the
#    discrete points.
# 
#    Nonlinear programming formulation
#        min   u     s.t.  u - phi >= 0, u + phi >= 0
#        x,u
# 
#    Specific problem: I = [-0.5,0.5]
#    phi(x,w) = 1/1+w - x1 exp(w x3) - x2 exp(w x4)
# 
#    Source: K. Oettershagen
#    "Ein superlinear knonvergenter algorithmus zur losung 
#     semi-infiniter optimierungsproblem",
#     Ph.D thesis, Bonn University, 1982
# 
#    SIF input: Nick Gould, February, 1994.
# 
#    classification = "C-CLOR2-AN-5-V"
# 
#    Discretization
# 
# IE M                   2
# IE M                   100
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "OET6"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling OET6.")
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
        v_["M"] = 500
        v_["LOWER"] = -0.5
        v_["UPPER"] = 0.5
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["DIFF"] = v_["UPPER"]-v_["LOWER"]
        v_["RM"] = Float64(v_["M"])
        v_["H"] = v_["DIFF"]/v_["RM"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["U"])
        push!(valA,Float64(1.0))
        for I = Int64(v_["0"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("LO"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"LO"*string(I))
            push!(irA,ig)
            push!(icA,ix_["U"])
            push!(valA,Float64(1.0))
            ig,ig_,_ = s2mpj_ii("UP"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"UP"*string(I))
            push!(irA,ig)
            push!(icA,ix_["U"])
            push!(valA,Float64(1.0))
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
        for I = Int64(v_["0"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["W"] = v_["RI"]*v_["H"]
            v_["W"] = v_["W"]+v_["LOWER"]
            v_["1+W"] = 1.0+v_["W"]
            v_["1/1+W"] = 1.0/v_["1+W"]
            v_["-1/1+W"] = -1.0*v_["1/1+W"]
            pbm.gconst[ig_["LO"*string(I)]] = Float64(v_["-1/1+W"])
            pbm.gconst[ig_["UP"*string(I)]] = Float64(v_["1/1+W"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXEYW", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["W"] = v_["RI"]*v_["H"]
            v_["W"] = v_["W"]+v_["LOWER"]
            ename = "E"*string(Int64(v_["1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXEYW")
            arrset(ielftype,ie,iet_["eXEYW"])
            ename = "E"*string(Int64(v_["1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["1"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="W",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["W"]))
            ename = "E"*string(Int64(v_["2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXEYW")
            arrset(ielftype,ie,iet_["eXEYW"])
            ename = "E"*string(Int64(v_["2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "E"*string(Int64(v_["2"]))*","*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            posep = findfirst(x->x=="W",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["W"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["M"])
            ig = ig_["LO"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["UP"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["1"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(Int64(v_["2"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLOR2-AN-5-V"
        pb.x0          = zeros(Float64,pb.n)
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eXEYW"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EYW = exp(EV_[2]*pbm.elpar[iel_][1])
        f_   = EV_[1]*EYW
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EYW
            g_[2] = EV_[1]*pbm.elpar[iel_][1]*EYW
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]*EYW
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]*pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EYW
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

