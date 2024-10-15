function OSBORNE2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OSBORNE2
#    *********
# 
#    Osborne second problem in 11 variables. This is a nonlinear equation version
#    of problem OSBORNEB.
# 
#    This function  is a nonlinear least squares with 65 groups.  Each
#    group has 4 nonlinear elements.
# 
#    Source:  Problem 19 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#32 (p.78).
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "C-NOR2-MN-11-65"
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "OSBORNE2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 65
        v_["N"] = 11
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(I))
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
        pbm.gconst[ig_["G1"]] = Float64(1.366)
        pbm.gconst[ig_["G2"]] = Float64(1.191)
        pbm.gconst[ig_["G3"]] = Float64(1.112)
        pbm.gconst[ig_["G4"]] = Float64(1.013)
        pbm.gconst[ig_["G5"]] = Float64(0.991)
        pbm.gconst[ig_["G6"]] = Float64(0.885)
        pbm.gconst[ig_["G7"]] = Float64(0.831)
        pbm.gconst[ig_["G8"]] = Float64(0.847)
        pbm.gconst[ig_["G9"]] = Float64(0.786)
        pbm.gconst[ig_["G10"]] = Float64(0.725)
        pbm.gconst[ig_["G11"]] = Float64(0.746)
        pbm.gconst[ig_["G12"]] = Float64(0.679)
        pbm.gconst[ig_["G13"]] = Float64(0.608)
        pbm.gconst[ig_["G14"]] = Float64(0.655)
        pbm.gconst[ig_["G15"]] = Float64(0.616)
        pbm.gconst[ig_["G16"]] = Float64(0.606)
        pbm.gconst[ig_["G17"]] = Float64(0.602)
        pbm.gconst[ig_["G18"]] = Float64(0.626)
        pbm.gconst[ig_["G19"]] = Float64(0.651)
        pbm.gconst[ig_["G20"]] = Float64(0.724)
        pbm.gconst[ig_["G21"]] = Float64(0.649)
        pbm.gconst[ig_["G22"]] = Float64(0.649)
        pbm.gconst[ig_["G23"]] = Float64(0.694)
        pbm.gconst[ig_["G24"]] = Float64(0.644)
        pbm.gconst[ig_["G25"]] = Float64(0.624)
        pbm.gconst[ig_["G26"]] = Float64(0.661)
        pbm.gconst[ig_["G27"]] = Float64(0.612)
        pbm.gconst[ig_["G28"]] = Float64(0.558)
        pbm.gconst[ig_["G29"]] = Float64(0.533)
        pbm.gconst[ig_["G30"]] = Float64(0.495)
        pbm.gconst[ig_["G31"]] = Float64(0.500)
        pbm.gconst[ig_["G32"]] = Float64(0.423)
        pbm.gconst[ig_["G33"]] = Float64(0.395)
        pbm.gconst[ig_["G34"]] = Float64(0.375)
        pbm.gconst[ig_["G35"]] = Float64(0.372)
        pbm.gconst[ig_["G36"]] = Float64(0.391)
        pbm.gconst[ig_["G37"]] = Float64(0.396)
        pbm.gconst[ig_["G38"]] = Float64(0.405)
        pbm.gconst[ig_["G39"]] = Float64(0.428)
        pbm.gconst[ig_["G40"]] = Float64(0.429)
        pbm.gconst[ig_["G41"]] = Float64(0.523)
        pbm.gconst[ig_["G42"]] = Float64(0.562)
        pbm.gconst[ig_["G43"]] = Float64(0.607)
        pbm.gconst[ig_["G44"]] = Float64(0.653)
        pbm.gconst[ig_["G45"]] = Float64(0.672)
        pbm.gconst[ig_["G46"]] = Float64(0.708)
        pbm.gconst[ig_["G47"]] = Float64(0.633)
        pbm.gconst[ig_["G48"]] = Float64(0.668)
        pbm.gconst[ig_["G49"]] = Float64(0.645)
        pbm.gconst[ig_["G50"]] = Float64(0.632)
        pbm.gconst[ig_["G51"]] = Float64(0.591)
        pbm.gconst[ig_["G52"]] = Float64(0.559)
        pbm.gconst[ig_["G53"]] = Float64(0.597)
        pbm.gconst[ig_["G54"]] = Float64(0.625)
        pbm.gconst[ig_["G55"]] = Float64(0.739)
        pbm.gconst[ig_["G56"]] = Float64(0.710)
        pbm.gconst[ig_["G57"]] = Float64(0.729)
        pbm.gconst[ig_["G58"]] = Float64(0.720)
        pbm.gconst[ig_["G59"]] = Float64(0.636)
        pbm.gconst[ig_["G60"]] = Float64(0.581)
        pbm.gconst[ig_["G61"]] = Float64(0.428)
        pbm.gconst[ig_["G62"]] = Float64(0.292)
        pbm.gconst[ig_["G63"]] = Float64(0.162)
        pbm.gconst[ig_["G64"]] = Float64(0.098)
        pbm.gconst[ig_["G65"]] = Float64(0.054)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(1.3)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(1.3)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(0.65)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(0.65)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(0.65)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(0.65)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(0.7)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(0.7)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(0.6)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(0.6)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(5.0)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(5.0)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(7.0)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(7.0)
        end
        if haskey(ix_,"X9")
            pb.x0[ix_["X9"]] = Float64(2.0)
        else
            pb.y0[findfirst(x->x==ig_["X9"],pbm.congrps)] = Float64(2.0)
        end
        if haskey(ix_,"X10")
            pb.x0[ix_["X10"]] = Float64(4.5)
        else
            pb.y0[findfirst(x->x==ig_["X10"],pbm.congrps)] = Float64(4.5)
        end
        if haskey(ix_,"X11")
            pb.x0[ix_["X11"]] = Float64(5.5)
        else
            pb.y0[findfirst(x->x==ig_["X11"],pbm.congrps)] = Float64(5.5)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePEXP", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        it,iet_,_ = s2mpj_ii( "ePEXP3", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftp,it,1,"T3")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["I-1"] = 1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["TI"] = 0.1*v_["RI-1"]
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePEXP")
            arrset(ielftype,ie,iet_["ePEXP"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePEXP3")
            arrset(ielftype,ie,iet_["ePEXP3"])
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X9"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePEXP3")
            arrset(ielftype,ie,iet_["ePEXP3"])
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X10"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
            ename = "D"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePEXP3")
            arrset(ielftype,ie,iet_["ePEXP3"])
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X11"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X8"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T3",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TI"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.04013774
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
        pb.pbclass = "C-NOR2-MN-11-65"
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

    elseif action == "ePEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPA = exp(-pbm.elpar[iel_][1]*EV_[2])
        FVAL = EV_[1]*EXPA
        f_   = FVAL
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPA
            g_[2] = -pbm.elpar[iel_][1]*FVAL
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -pbm.elpar[iel_][1]*EXPA
                H_[2,1] = H_[1,2]
                H_[2,2] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*FVAL
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePEXP3"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TMV2 = pbm.elpar[iel_][1]-EV_[2]
        TMV2SQ = TMV2*TMV2
        EXPA = exp(-TMV2SQ*EV_[3])
        FVAL = EV_[1]*EXPA
        A = 2.0*TMV2*EV_[3]
        f_   = FVAL
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPA
            g_[2] = A*FVAL
            g_[3] = -TMV2SQ*FVAL
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = A*EXPA
                H_[2,1] = H_[1,2]
                H_[1,3] = -TMV2SQ*EXPA
                H_[3,1] = H_[1,3]
                H_[2,2] = (A*A-2.0*EV_[3])*FVAL
                H_[2,3] = (2.0*TMV2-A*TMV2SQ)*FVAL
                H_[3,2] = H_[2,3]
                H_[3,3] = TMV2SQ*TMV2SQ*FVAL
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

