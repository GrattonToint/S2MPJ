function RAT43(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : RAT43
#    *********
# 
#    NIST Data fitting problem RAT43 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = b1 / ((1+exp[b2-b3*x])**(1/b4)) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Ratkowsky, D.A. (1983).  
#      Nonlinear Regression Modeling.
#      New York, NY:  Marcel Dekker, pp. 62 and 88.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-NOR2-MN-4-15"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "RAT43"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 15
        v_["N"] = 4
        v_["1"] = 1
        v_["X1"] = 1.0
        v_["X2"] = 2.0
        v_["X3"] = 3.0
        v_["X4"] = 4.0
        v_["X5"] = 5.0
        v_["X6"] = 6.0
        v_["X7"] = 7.0
        v_["X8"] = 8.0
        v_["X9"] = 9.0
        v_["X10"] = 10.0
        v_["X11"] = 11.0
        v_["X12"] = 12.0
        v_["X13"] = 13.0
        v_["X14"] = 14.0
        v_["X15"] = 15.0
        v_["Y1"] = 16.08
        v_["Y2"] = 33.83
        v_["Y3"] = 65.80
        v_["Y4"] = 97.20
        v_["Y5"] = 191.55
        v_["Y6"] = 326.20
        v_["Y7"] = 386.87
        v_["Y8"] = 520.53
        v_["Y9"] = 590.03
        v_["Y10"] = 651.92
        v_["Y11"] = 724.93
        v_["Y12"] = 699.56
        v_["Y13"] = 689.96
        v_["Y14"] = 637.56
        v_["Y15"] = 717.41
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
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
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(100.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(100.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(10.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(10.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"B4")
            pb.x0[ix_["B4"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["B4"],pbm.congrps)] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE")
            arrset(ielftype,ie,iet_["eE"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
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
        pb.pbclass = "C-NOR2-MN-4-15"
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

    elseif action == "eE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V2MV3X = EV_[2]-EV_[3]*pbm.elpar[iel_][1]
        V4INV = 1.0/EV_[4]
        V4INVP = V4INV+1.0
        E = exp(V2MV3X)
        E2 = E*E
        EP1 = E+1.0
        EP1L = log(EP1)
        EP14 = EP1^V4INV
        EP14P1 = EP1^V4INVP
        EP14P2 = EP1^(V4INV+2.0)
        VE = EV_[4]*EP14P1
        VE2 = EV_[4]*EP14P2
        V42EPP = EP14*EV_[4]^2
        V42EP2 = EP14P1*EV_[4]^2
        V42EP3 = EP14P1*EV_[4]^3
        f_   = EV_[1]/EP14
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EP14
            g_[2] = -EV_[1]*E/VE
            g_[3] = EV_[1]*pbm.elpar[iel_][1]*E/VE
            g_[4] = EV_[1]*EP1L/V42EPP
            if nargout>2
                H_ = zeros(Float64,4,4)
                H_[1,2] = -E/VE
                H_[2,1] = H_[1,2]
                H_[1,3] = pbm.elpar[iel_][1]*E/VE
                H_[3,1] = H_[1,3]
                H_[1,4] = EP1L/V42EPP
                H_[4,1] = H_[1,4]
                H_[2,2] = EV_[1]*(E2*V4INVP/VE2-E/VE)
                H_[2,3] = EV_[1]*pbm.elpar[iel_][1]*(E/VE-E2*V4INVP/VE2)
                H_[3,2] = H_[2,3]
                H_[2,4] = EV_[1]*E*(1.0/V42EP2-EP1L/V42EP3)
                H_[4,2] = H_[2,4]
                H_[3,3] = EV_[1]*pbm.elpar[iel_][1]^2*(E2*V4INVP/VE2-E/VE)
                H_[3,4] = EV_[1]*pbm.elpar[iel_][1]*E*(EP1L/V42EP3-1.0/V42EP2)
                H_[4,3] = H_[3,4]
                H_[4,4] = (EV_[1]/EP14)*(EP1L^2/EV_[4]^4-2.0*EP1L/EV_[4]^3)
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

