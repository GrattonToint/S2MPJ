function YFITNE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This problem arises in measuring angles and distances to a vibrating beam
#    using a laser-Doppler velocimeter.
#    This is a nonlinear equation variant of the bounded constrained
#    problem YFIT.
# 
#    Source:
#    an exercize for L. Watson course on LANCELOT in the Spring 1993.
# 
#    SIF input: B. E. Lindholm, Virginia Tech., Spring 1993,
#               modified by Ph. Toint, March 1994.
#               derivatives corrected by Nick Gould, June 2019.
# 
#    classification = "C-NOR2-MN-3-17"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "YFITNE"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["zero"] = 0
        v_["p"] = 16
        v_["realp"] = 16.0
        v_["y0"] = 21.158931
        v_["y1"] = 17.591719
        v_["y2"] = 14.046854
        v_["y3"] = 10.519732
        v_["y4"] = 7.0058392
        v_["y5"] = 3.5007293
        v_["y6"] = 0.0000000
        v_["y7"] = -3.5007293
        v_["y8"] = -7.0058392
        v_["y9"] = -10.519732
        v_["y10"] = -14.046854
        v_["y11"] = -17.591719
        v_["y12"] = -21.158931
        v_["y13"] = -24.753206
        v_["y14"] = -28.379405
        v_["y15"] = -32.042552
        v_["y16"] = -35.747869
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("alpha",ix_)
        arrset(pb.xnames,iv,"alpha")
        iv,ix_,_ = s2mpj_ii("beta",ix_)
        arrset(pb.xnames,iv,"beta")
        iv,ix_,_ = s2mpj_ii("dist",ix_)
        arrset(pb.xnames,iv,"dist")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for i = Int64(v_["zero"]):Int64(v_["p"])
            ig,ig_,_ = s2mpj_ii("diff"*string(i),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"diff"*string(i))
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
        for i = Int64(v_["zero"]):Int64(v_["p"])
            pbm.gconst[ig_["diff"*string(i)]] = Float64(v_["y"*string(i)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["alpha"]] = Float64(0.60)
        pb.x0[ix_["beta"]] = Float64(-0.60)
        pb.x0[ix_["dist"]] = Float64(20.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "etanab", iet_)
        loaset(elftv,it,1,"a1")
        loaset(elftv,it,2,"b1")
        loaset(elftv,it,3,"d1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"point")
        loaset(elftp,it,2,"count")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for i = Int64(v_["zero"]):Int64(v_["p"])
            v_["index"] = Float64(i)
            ename = "est"*string(i)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"etanab")
            arrset(ielftype,ie,iet_["etanab"])
            vname = "alpha"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="a1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "beta"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="b1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "dist"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="d1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="point",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["index"]))
            posep = findfirst(x->x=="count",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["realp"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for i = Int64(v_["zero"]):Int64(v_["p"])
            ig = ig_["diff"*string(i)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["est"*string(i)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
# LO SOLUTION            0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-MN-3-17"
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

    elseif action == "etanab"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        frac = pbm.elpar[iel_][1]/pbm.elpar[iel_][2]
        ttan = tan(EV_[1]*(1.0-frac)+EV_[2]*frac)
        tsec = 1.0/cos(EV_[1]*(1.0-frac)+EV_[2]*frac)
        tsec2 = tsec*tsec
        f_   = EV_[3]*ttan
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[3]*(1.0-frac)*tsec2
            g_[2] = EV_[3]*frac*tsec2
            g_[3] = ttan
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0*EV_[3]*((1.0-frac)^2)*tsec2*ttan
                H_[2,2] = 2.0*EV_[3]*(frac^2)*tsec2*ttan
                H_[1,2] = 2.0*EV_[3]*(1.0-frac)*frac*tsec2*ttan
                H_[2,1] = H_[1,2]
                H_[1,3] = (1.0-frac)*tsec2
                H_[3,1] = H_[1,3]
                H_[2,3] = frac*tsec2
                H_[3,2] = H_[2,3]
                H_[3,3] = 0.0
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

