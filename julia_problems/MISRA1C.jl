function MISRA1C(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : MISRA1C
#    *********
# 
#    NIST Data fitting problem MISRA1C given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 * (1-(1+2*b2*x)**(-.5)) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Misra, D., NIST (1978).  
#      Dental Research Monomolecular Adsorption Study.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-CNOR2-MN-2-14"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "MISRA1C"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling MISRA1C.")
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
        v_["M"] = 14
        v_["N"] = 2
        v_["1"] = 1
        v_["X1"] = 77.6
        v_["X2"] = 114.9
        v_["X3"] = 141.1
        v_["X4"] = 190.8
        v_["X5"] = 239.9
        v_["X6"] = 289.0
        v_["X7"] = 332.8
        v_["X8"] = 378.4
        v_["X9"] = 434.8
        v_["X10"] = 477.3
        v_["X11"] = 536.8
        v_["X12"] = 593.1
        v_["X13"] = 689.1
        v_["X14"] = 760.0
        v_["Y1"] = 10.07
        v_["Y2"] = 14.73
        v_["Y3"] = 17.94
        v_["Y4"] = 23.93
        v_["Y5"] = 29.61
        v_["Y6"] = 35.18
        v_["Y7"] = 40.02
        v_["Y8"] = 44.82
        v_["Y9"] = 50.76
        v_["Y10"] = 55.05
        v_["Y11"] = 61.01
        v_["Y12"] = 66.40
        v_["Y13"] = 75.47
        v_["Y14"] = 81.78
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
            push!(irA,ig)
            push!(icA,ix_["B1"])
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(500.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(500.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(0.0001)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(0.0001)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE3", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE3")
            arrset(ielftype,ie,iet_["eE3"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CNOR2-MN-2-14"
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

    elseif action == "eE3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        R = 1.0+2.0*EV_[2]*pbm.elpar[iel_][1]
        ROOTR = sqrt(R)
        R3 = ROOTR^3
        f_   = -EV_[1]/ROOTR
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -1.0/ROOTR
            g_[2] = EV_[1]*pbm.elpar[iel_][1]/R3
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = pbm.elpar[iel_][1]/R3
                H_[2,1] = H_[1,2]
                H_[2,2] = -3.0*EV_[1]*pbm.elpar[iel_][1]^2/ROOTR^5
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

