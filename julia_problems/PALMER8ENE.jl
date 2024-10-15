function PALMER8ENE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER8ENE
#    *********
# 
#    A nonlinear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + L * EXP( -K X**2 )
# 
#    Source:
#    M.  Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "C-NOR2-RN-8-12"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER8ENE"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 23
        v_["1"] = 1
        v_["12"] = 12
        v_["X12"] = 0.000000
        v_["X13"] = 0.174533
        v_["X14"] = 0.314159
        v_["X15"] = 0.436332
        v_["X16"] = 0.514504
        v_["X17"] = 0.610865
        v_["X18"] = 0.785398
        v_["X19"] = 0.959931
        v_["X20"] = 1.134464
        v_["X21"] = 1.308997
        v_["X22"] = 1.483530
        v_["X23"] = 1.570796
        v_["Y12"] = 4.757534
        v_["Y13"] = 3.121416
        v_["Y14"] = 1.207606
        v_["Y15"] = 0.131916
        v_["Y16"] = 0.000000
        v_["Y17"] = 0.258514
        v_["Y18"] = 3.380161
        v_["Y19"] = 10.762813
        v_["Y20"] = 23.745996
        v_["Y21"] = 44.471864
        v_["Y22"] = 76.541947
        v_["Y23"] = 97.874528
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("A0",ix_)
        arrset(pb.xnames,iv,"A0")
        iv,ix_,_ = s2mpj_ii("A2",ix_)
        arrset(pb.xnames,iv,"A2")
        iv,ix_,_ = s2mpj_ii("A4",ix_)
        arrset(pb.xnames,iv,"A4")
        iv,ix_,_ = s2mpj_ii("A6",ix_)
        arrset(pb.xnames,iv,"A6")
        iv,ix_,_ = s2mpj_ii("A8",ix_)
        arrset(pb.xnames,iv,"A8")
        iv,ix_,_ = s2mpj_ii("A10",ix_)
        arrset(pb.xnames,iv,"A10")
        iv,ix_,_ = s2mpj_ii("K",ix_)
        arrset(pb.xnames,iv,"K")
        iv,ix_,_ = s2mpj_ii("L",ix_)
        arrset(pb.xnames,iv,"L")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["12"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            v_["XQUART"] = v_["XSQR"]*v_["XSQR"]
            v_["X**6"] = v_["XSQR"]*v_["XQUART"]
            v_["X**8"] = v_["XSQR"]*v_["X**6"]
            v_["X**10"] = v_["XSQR"]*v_["X**8"]
            v_["X**12"] = v_["XSQR"]*v_["X**10"]
            v_["X**14"] = v_["XSQR"]*v_["X**12"]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"O"*string(I))
            iv = ix_["A0"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["A2"]
            pbm.A[ig,iv] += Float64(v_["XSQR"])
            iv = ix_["A4"]
            pbm.A[ig,iv] += Float64(v_["XQUART"])
            iv = ix_["A6"]
            pbm.A[ig,iv] += Float64(v_["X**6"])
            iv = ix_["A8"]
            pbm.A[ig,iv] += Float64(v_["X**8"])
            iv = ix_["A10"]
            pbm.A[ig,iv] += Float64(v_["X**10"])
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
        for I = Int64(v_["12"]):Int64(v_["M"])
            pbm.gconst[ig_["O"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["A0"]] = -Inf
        pb.xupper[ix_["A0"]] = +Inf
        pb.xlower[ix_["A2"]] = -Inf
        pb.xupper[ix_["A2"]] = +Inf
        pb.xlower[ix_["A4"]] = -Inf
        pb.xupper[ix_["A4"]] = +Inf
        pb.xlower[ix_["A6"]] = -Inf
        pb.xupper[ix_["A6"]] = +Inf
        pb.xlower[ix_["A8"]] = -Inf
        pb.xupper[ix_["A8"]] = +Inf
        pb.xlower[ix_["A10"]] = -Inf
        pb.xupper[ix_["A10"]] = +Inf
        pb.xlower[ix_["L"]] = -Inf
        pb.xupper[ix_["L"]] = +Inf
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"K")
        loaset(elftv,it,2,"L")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"XSQR")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["12"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "K"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="K",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "L"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="L",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="XSQR",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["XSQR"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["12"]):Int64(v_["M"])
            ig = ig_["O"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
# LO PALMER8E                0.0
#    Solution
# LO SOLTN              1.48003482D-04
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
        pb.pbclass = "C-NOR2-RN-8-12"
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

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPON = exp(-EV_[1]*pbm.elpar[iel_][1])
        f_   = EV_[2]*EXPON
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.elpar[iel_][1]*EV_[2]*EXPON
            g_[2] = EXPON
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EV_[2]*EXPON
                H_[1,2] = -pbm.elpar[iel_][1]*EXPON
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

