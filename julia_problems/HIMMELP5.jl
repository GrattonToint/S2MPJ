function HIMMELP5(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HIMMELP5
#    *********
# 
#    A nonlinear problem with inequality constraints, attributed to Himmelblau
#    by B.N. Pshenichnyj (case IV).
# 
#    The problem is nonconvex and has redundant constraints at the solution.
# 
#    Source: 
#    B.N. Pshenichnyj
#    "The Linearization Method for Constrained Optimization",
#    Springer Verlag, SCM Series 22, Heidelberg, 1994
# 
#    SIF input: Ph. Toint, December 1994.
# 
#    classification = "C-OQR2-AN-2-3"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HIMMELP5"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["B1"] = 0.1963666677
        v_["B1"] = 75.0+v_["B1"]
        v_["B2"] = -.8112755343
        v_["B2"] = -3.0+v_["B2"]
        v_["B6"] = -.8306567613
        v_["B6"] = -6.0+v_["B6"]
        v_["-B2"] = -1.0*v_["B2"]
        v_["-B6"] = -1.0*v_["B6"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(v_["-B2"])
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(v_["-B6"])
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C5")
        ig,ig_,_ = s2mpj_ii("C8",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C8")
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("C9",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C9")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(100.0)
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
        pbm.gconst[ig_["OBJ"]] = Float64(v_["B1"])
        pbm.gconst[ig_["C5"]] = Float64(700.0)
        pbm.gconst[ig_["C9"]] = Float64(2775.0)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 54.0
        pb.xupper[ix_["X1"]] = 75.0
        pb.xupper[ix_["X2"]] = 65.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X1"]] = Float64(68.8)
        pb.x0[ix_["X2"]] = Float64(31.2)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eOBNL", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OB"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOBNL")
        arrset(ielftype,ie,iet_["eOBNL"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X1X2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en2PR")
        arrset(ielftype,ie,iet_["en2PR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X1SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X2SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OB"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X1X2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C8"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X1SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.008))
        ig = ig_["C9"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X2SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN                -59.01312394
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = fill(Inf,pb.nge)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OQR2-AN-2-3"
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

    elseif action == "eOBNL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B3 = .1269366345
        B4 = -0.20567665
        B4 = 0.01*B4
        B5 = 0.103450e-4
        B7 = .0302344793
        B8 = -0.12813448
        B8 = 0.01*B8
        B9 = 0.352599e-4
        B10 = -0.2266e-6
        B11 = 0.2564581253
        B12 = -.003460403
        B13 = 0.135139e-4
        B14 = -.1064434908
        B14 = B14-28.0
        B15 = -0.52375e-5
        B16 = -0.63e-8
        B17 = 0.7e-9
        B18 = 0.3405462
        B18 = 0.001*B18
        B19 = -0.16638e-5
        B20 = -2.86731123
        B20 = B20-0.92e-8
        A = B7*EV_[1]+B8*EV_[1]^2+B9*EV_[1]^3+B10*EV_[1]^4
        DADX = B7+2.0*B8*EV_[1]+3.0*B9*EV_[1]^2+4.0*B10*EV_[1]^3
        D2ADXX = 2.0*B8+6.0*B9*EV_[1]+12.0*B10*EV_[1]^2
        B = B18*EV_[1]+B15*EV_[1]^2+B16*EV_[1]^3
        DBDX = B18+2.0*B15*EV_[1]+3.0*B16*EV_[1]^2
        D2BDXX = 2.0*B15+6.0*B16*EV_[1]
        C = B3*EV_[1]^2+B4*EV_[1]^3+B5*EV_[1]^4
        DCDX = 2.0*B3*EV_[1]+3.0*B4*EV_[1]^2+4.0*B5*EV_[1]^3
        D2CDXX = 2.0*B3+6.0*B4*EV_[1]+12.0*B5*EV_[1]^2
        F = B11*EV_[2]^2+B12*EV_[2]^3+B13*EV_[2]^4
        DFDY = 2.0*B11*EV_[2]+3.0*B12*EV_[2]^2+4.0*B13*EV_[2]^3
        D2FDYY = 2.0*B11+6.0*B12*EV_[2]+12.0*B13*EV_[2]^2
        G = B17*EV_[1]^3+B19*EV_[1]
        DGDX = B19+3.0*B17*EV_[1]^2
        D2GDXX = 6.0*B17*EV_[1]
        E = exp(0.0005*EV_[1]*EV_[2])
        DEDX = 0.0005*EV_[2]*E
        DEDY = 0.0005*EV_[1]*E
        D2EDXX = 0.0005*EV_[2]*DEDX
        D2EDXY = 0.0005*(EV_[2]*DEDY+E)
        D2EDYY = 0.0005*EV_[1]*DEDY
        f_   = C+EV_[2]*A+F+B14/(1.0+EV_[2])+B*EV_[2]^2+G*EV_[2]^3+B20*E
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DCDX+EV_[2]*DADX+DBDX*EV_[2]^2+DGDX*EV_[2]^3+B20*DEDX
            g_[2] = A+DFDY-B14/(1.0+EV_[2])^2+2.0*B*EV_[2]+3.0*G*EV_[2]^2+B20*DEDY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = D2CDXX+EV_[2]*D2ADXX+D2BDXX*EV_[2]^2+D2GDXX*EV_[2]^3+B20*D2EDXX
                H_[1,2] = DADX+2.0*EV_[2]*DBDX+3.0*DGDX*EV_[2]^2+B20*D2EDXY
                H_[2,1] = H_[1,2]
                H_[2,2] = D2FDYY+2.0*B14/(1.0+EV_[2])^3+2.0*B+6.0*G*EV_[2]+B20*D2EDYY
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "en2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
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

