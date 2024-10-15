function PALMER3(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER3
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#     model: H-N=C=S TZVP + MP2
#    fitting Y to A X**2 + B / ( C + X**2 / D ), B, C, D nonnegative.
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "C-SBR2-RN-4-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER3"

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
        v_["X1"] = -1.658063
        v_["X2"] = -1.570796
        v_["X3"] = -1.396263
        v_["X4"] = -1.221730
        v_["X5"] = -1.047198
        v_["X6"] = -0.872665
        v_["X7"] = -0.766531
        v_["X8"] = -0.698132
        v_["X9"] = -0.523599
        v_["X10"] = -0.349066
        v_["X11"] = -0.174533
        v_["X12"] = 0.0
        v_["X13"] = 0.174533
        v_["X14"] = 0.349066
        v_["X15"] = 0.523599
        v_["X16"] = 0.698132
        v_["X17"] = 0.766531
        v_["X18"] = 0.872665
        v_["X19"] = 1.047198
        v_["X20"] = 1.221730
        v_["X21"] = 1.396263
        v_["X22"] = 1.570796
        v_["X23"] = 1.658063
        v_["Y1"] = 64.87939
        v_["Y2"] = 50.46046
        v_["Y3"] = 28.2034
        v_["Y4"] = 13.4575
        v_["Y5"] = 4.6547
        v_["Y6"] = 0.59447
        v_["Y7"] = 0.0000
        v_["Y8"] = 0.2177
        v_["Y9"] = 2.3029
        v_["Y10"] = 5.5191
        v_["Y11"] = 8.5519
        v_["Y12"] = 9.8919
        v_["Y13"] = 8.5519
        v_["Y14"] = 5.5191
        v_["Y15"] = 2.3029
        v_["Y16"] = 0.2177
        v_["Y17"] = 0.0000
        v_["Y18"] = 0.59447
        v_["Y19"] = 4.6547
        v_["Y20"] = 13.4575
        v_["Y21"] = 28.2034
        v_["Y22"] = 50.46046
        v_["Y23"] = 64.87939
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("A",ix_)
        arrset(pb.xnames,iv,"A")
        iv,ix_,_ = s2mpj_ii("B",ix_)
        arrset(pb.xnames,iv,"B")
        iv,ix_,_ = s2mpj_ii("C",ix_)
        arrset(pb.xnames,iv,"C")
        iv,ix_,_ = s2mpj_ii("D",ix_)
        arrset(pb.xnames,iv,"D")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["A"]
            pbm.A[ig,iv] += Float64(v_["XSQR"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["O"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["A"]] = -Inf
        pb.xupper[ix_["A"]] = +Inf
        pb.xlower[ix_["B"]] = 0.00001
        pb.xlower[ix_["C"]] = 0.00001
        pb.xlower[ix_["D"]] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eQUOT", iet_)
        loaset(elftv,it,1,"B")
        loaset(elftv,it,2,"C")
        loaset(elftv,it,3,"D")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"XSQR")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eQUOT")
            arrset(ielftype,ie,iet_["eQUOT"])
            vname = "B"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="B",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="C",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "D"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="D",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="XSQR",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["XSQR"]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["O"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
#    SOLTN               2265.95822
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-RN-4-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eQUOT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DENOM = 1.0/(EV_[2]+pbm.elpar[iel_][1]/EV_[3])
        f_   = EV_[1]*DENOM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DENOM
            g_[2] = -EV_[1]*DENOM*DENOM
            g_[3] = EV_[1]*pbm.elpar[iel_][1]*(DENOM/EV_[3])^2
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = -DENOM*DENOM
                H_[2,1] = H_[1,2]
                H_[1,3] = pbm.elpar[iel_][1]*(DENOM/EV_[3])^2
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*EV_[1]*DENOM^3
                H_[2,3] = -2.0*EV_[1]*pbm.elpar[iel_][1]*DENOM*(DENOM/EV_[3])^2
                H_[3,2] = H_[2,3]
                H_[3,3]  = (
                      2.0*EV_[1]*DENOM^3*(pbm.elpar[iel_][1]/EV_[3]^2)^2-2.0*EV_[1]*DENOM^2*pbm.elpar[iel_][1]/EV_[3]^3)
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

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

