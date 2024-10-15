function HS54(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS54
#    *********
# 
#    Source: problem 54, incorrectly stated in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    Betts problem 11.7, JOTA 21, 1977, pp.137-174.
#    SIF input: A.R. Conn, April 1990 and Nick Gould, October 1990
# 
#    classification = "C-OLR2-AN-6-1"
# 
#    some useful parameters, including N, the number of variables.
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS54"

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
        v_["1"] = 1
        v_["6"] = 6
        v_["RHO"] = 2.0e-1
        v_["RHOSQR"] = v_["RHO"]*v_["RHO"]
        v_["1-RHOSQR"] = 1.0-v_["RHOSQR"]
        v_["FACTOR"] = 1.0/v_["1-RHOSQR"]
        v_["MU1"] = 1.0e+4
        v_["MU2"] = 1.0e+0
        v_["MU3"] = 2.0e+6
        v_["MU4"] = 1.0e+1
        v_["MU5"] = 1.0e-3
        v_["MU6"] = 1.0e+8
        v_["SIGMA1"] = 8.0e+3
        v_["SIGMA2"] = 1.0e+0
        v_["SIGMA3"] = 7.0e+6
        v_["SIGMA4"] = 5.0e+1
        v_["SIGMA5"] = 5.0e-2
        v_["SIGMA6"] = 5.0e+8
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
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CON1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"CON1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(4.0e+3)
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
        v_["0.2SI1"] = 2.0e-1*v_["SIGMA1"]
        v_["2000SI2"] = 2.0e+3*v_["SIGMA2"]
        v_["4000MU2"] = 4.0e+3*v_["MU2"]
        v_["RHS"] = v_["MU1"]+v_["4000MU2"]
        v_["RHS"] = v_["RHS"]+v_["0.2SI1"]
        v_["RHS"] = v_["RHS"]+v_["2000SI2"]
        pbm.gconst[ig_["CON1"]] = Float64(v_["RHS"])
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["X1"]] = 2.0e+4
        pb.xlower[ix_["X2"]] = - 1.0e+1
        pb.xupper[ix_["X2"]] = 1.0e+1
        pb.xupper[ix_["X3"]] = 1.0e+7
        pb.xupper[ix_["X4"]] = 2.0e+1
        pb.xlower[ix_["X5"]] = - 1.0e+0
        pb.xupper[ix_["X5"]] = 1.0e+0
        pb.xupper[ix_["X6"]] = 2.0e+8
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(6.0e+3)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(6.0e+3)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(1.5e+0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(1.5e+0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(4.0e+6)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(4.0e+6)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(2.0e+0)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(2.0e+0)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(3.0e-3)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(3.0e-3)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(5.0e+7)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(5.0e+7)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"V1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"MU")
        loaset(elftp,it,2,"SIGMA")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftp,it,1,"MU1")
        loaset(elftp,it,2,"MU2")
        loaset(elftp,it,3,"SIGMA1")
        loaset(elftp,it,4,"SIGMA2")
        loaset(elftp,it,5,"RHO")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["6"])
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eSQR")
                arrset(ielftype,ie,iet_["eSQR"])
            end
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="MU",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["MU"*string(I)]))
            posep = findfirst(x->x=="SIGMA",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["SIGMA"*string(I)]))
        end
        ename = "F1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="RHO",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["RHO"]))
        posep = findfirst(x->x=="MU1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["MU1"]))
        posep = findfirst(x->x=="MU2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["MU2"]))
        posep = findfirst(x->x=="SIGMA1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["SIGMA1"]))
        posep = findfirst(x->x=="SIGMA2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["SIGMA2"]))
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gNORMAL",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        arrset(pbm.grftype,ig,"gNORMAL")
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["FACTOR"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["FACTOR"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0e+0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(1.0e+0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0e+0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E6"])
        loaset(pbm.grelw,ig,posel,Float64(1.0e+0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["F1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["FACTOR"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.90807482
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
        pb.pbclass = "C-OLR2-AN-6-1"
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

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V1MP = (EV_[1]-pbm.elpar[iel_][1])/pbm.elpar[iel_][2]
        f_   = V1MP^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*V1MP/pbm.elpar[iel_][2]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0/pbm.elpar[iel_][2]^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        TERM1 = (EV_[1]-pbm.elpar[iel_][1])/pbm.elpar[iel_][3]
        TERM2 = (EV_[2]-pbm.elpar[iel_][2])/pbm.elpar[iel_][4]
        RHO2 = pbm.elpar[iel_][5]+pbm.elpar[iel_][5]
        f_   = RHO2*TERM1*TERM2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = RHO2*TERM2/pbm.elpar[iel_][3]
            g_[2] = RHO2*TERM1/pbm.elpar[iel_][4]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = RHO2/(pbm.elpar[iel_][3]*pbm.elpar[iel_][4])
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gNORMAL"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPHV = exp(-0.5*GVAR_)
        f_= -EXPHV
        if nargout>1
            g_ = 5.0e-1*EXPHV
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = -2.5e-1*EXPHV
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

