using SpecialFunctions
function HS69(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS69
#    *********
# 
#    This is a cost optimal inspection plan.
# 
#    Source: problem 69 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-COOR2-MN-4-2"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 25 XI 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS69"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 4
        v_["A"] = 0.1
        v_["B"] = 1000.0
        v_["D"] = 1.0
        v_["NN"] = 4.0
        v_["1"] = 1
        v_["AN"] = v_["A"]*v_["NN"]
        v_["ROOTN"] = sqrt(v_["NN"])
        v_["DROOTN"] = v_["D"]*v_["ROOTN"]
        v_["-DROOTN"] = -v_["DROOTN"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C1")
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(1.0e+0))
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C2")
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(1.0e+0))
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
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 0.0001
        pb.xupper[ix_["X1"]] = 100.0
        pb.xupper[ix_["X2"]] = 100.0
        pb.xupper[ix_["X3"]] = 2.0
        pb.xupper[ix_["X4"]] = 2.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.y0 = fill(Float64(1.0),pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eRECIP", iet_)
        loaset(elftv,it,1,"X1")
        it,iet_,_ = s2mpj_ii( "eNASTYEXP", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X3")
        loaset(elftv,it,3,"X4")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"B")
        it,iet_,_ = s2mpj_ii( "ePHI", iet_)
        loaset(elftv,it,1,"X2")
        loaset(elftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eRECIP")
        arrset(ielftype,ie,iet_["eRECIP"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "OE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eNASTYEXP")
        arrset(ielftype,ie,iet_["eNASTYEXP"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="B",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["B"]))
        ename = "C1E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePHI")
        arrset(ielftype,ie,iet_["ePHI"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0e+0))
        ename = "C2E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePHI")
        arrset(ielftype,ie,iet_["ePHI"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["DROOTN"]))
        ename = "C2E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePHI")
        arrset(ielftype,ie,iet_["ePHI"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["-DROOTN"]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["AN"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0e+0))
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-2.0e+0))
        ig = ig_["C2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0e+0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C2E2"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0e+0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -9.20389D-01
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
        pb.pbclass = "C-COOR2-MN-4-2"
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

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,3.9894228040143270e-01)
        return pbm

    elseif action == "eRECIP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = 1.0e+0/EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -1.0e+0/EV_[1]^2
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0e+0/EV_[1]^3
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eNASTYEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E = exp(EV_[1])
        R = pbm.elpar[iel_][1]*(E-1.0e+0)-EV_[2]
        S = E-1.0e+0+EV_[3]
        F1 = -(-(R*EV_[3]/S))/(EV_[1]*EV_[1])
        F2 = -R/(S*EV_[1])
        F3 = -EV_[3]/(S*EV_[1])
        F4 = (EV_[3]*R)/(EV_[1]*S*S)
        F11 = 2.0e+0*(-(R*EV_[3]/S))/(EV_[1]*EV_[1]*EV_[1])
        F12 = E*EV_[3]*((R/S)-pbm.elpar[iel_][1])/(EV_[1]*S)
        F11 = F11+F12
        F12 = R/(S*EV_[1]*EV_[1])
        F13 = EV_[3]/(S*EV_[1]*EV_[1])
        F14 = R/(S*S*EV_[1])
        F15 = -1.0e+0/(S*EV_[1])
        F16 = -EV_[3]*R/(S*S*EV_[1]*EV_[1])
        F17 = EV_[3]/(S*S*EV_[1])
        F18 = -2.0e+0*EV_[3]*R/(EV_[1]*S*S*S)
        f_   = EV_[3]*R/(S*EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -F1-(pbm.elpar[iel_][1]*E*F3)-(E*F4)
            g_[2] = F3
            g_[3] = -F2-F4
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1]  = (
                      -F11-(2.0e+0*pbm.elpar[iel_][1]*E*F13)-(2.0e+0*E*F16)-(2.0e+0*pbm.elpar[iel_][1]*E*E*F17)-(E*E*F18))
                H_[2,1] = F13+(E*F17)
                H_[1,2] = H_[2,1]
                H_[3,1] = (-(E*F14)-(pbm.elpar[iel_][1]*E*F15)-F16-(pbm.elpar[iel_][1]*E*F17)-
                     (E*F18)-F12)
                H_[1,3] = H_[3,1]
                H_[3,2] = F15+F17
                H_[2,3] = H_[3,2]
                H_[3,3] = -(2.0e+0*F14)-F18
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePHI"

        x2      = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        mx2pp   = -x2+pbm.elpar[iel_][1]
        E       = exp(-0.5*mx2pp^2)
        thearg  =  7.071067811865475e-1 * abs( mx2pp )
        if ( mx2pp >= 0 )
            f_ = 0.5+0.5*erf( thearg )
        else
            f_ = 0.5-0.5*erf( thearg )
        end
        if nargout>1
            dim = length(x2)
            g_  = zeros(Float64,dim)
            g_[1] = -pbm.efpar[1]*E
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -pbm.efpar[1]*E*mx2pp
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv",
                       "cJxv","cJtxv","cIJtxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy",
                       "LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [1,0]
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

