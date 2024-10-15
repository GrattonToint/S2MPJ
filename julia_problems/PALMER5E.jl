function PALMER5E(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER5E
#    *********
# 
#    A nonlinear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 T_0 + A2 T_2 + A4 T_4 + A6 T_6 + A8 T_8 +
#                 A10 T_10 + A12 T_12 + A14 T_14
#                 + L * EXP( -K X**2 )
#    where T_i is the i-th (shifted) Chebyshev polynomial
# 
#    Source:
#    M.  Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
# 
#    classification = "C-SBR2-RN-8-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER5E"

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
        v_["2"] = 2
        v_["12"] = 12
        v_["14"] = 14
        v_["X12"] = 0.000000
        v_["X13"] = 1.570796
        v_["X14"] = 1.396263
        v_["X15"] = 1.308997
        v_["X16"] = 1.221730
        v_["X17"] = 1.125835
        v_["X18"] = 1.047198
        v_["X19"] = 0.872665
        v_["X20"] = 0.698132
        v_["X21"] = 0.523599
        v_["X22"] = 0.349066
        v_["X23"] = 0.174533
        v_["B"] = v_["X13"]
        v_["A"] = -1.0e+0*v_["B"]
        v_["DIFF"] = 2.0e+0*v_["B"]
        v_["Y12"] = 83.57418
        v_["Y13"] = 81.007654
        v_["Y14"] = 18.983286
        v_["Y15"] = 8.051067
        v_["Y16"] = 2.044762
        v_["Y17"] = 0.000000
        v_["Y18"] = 1.170451
        v_["Y19"] = 10.479881
        v_["Y20"] = 25.785001
        v_["Y21"] = 44.126844
        v_["Y22"] = 62.822177
        v_["Y23"] = 77.719674
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
            v_["T0"] = 1.0e+0
            v_["Y"] = 2.0e+0*v_["X"*string(I)]
            v_["Y"] = v_["Y"]-v_["A"]
            v_["Y"] = v_["Y"]-v_["B"]
            v_["Y"] = v_["Y"]/v_["DIFF"]
            v_["T1"] = v_["Y"]
            v_["2Y"] = 2.0e+0*v_["Y"]
            for J = Int64(v_["2"]):Int64(v_["14"])
                v_["J-1"] = -1+J
                v_["J-2"] = -2+J
                v_["T"*string(J)] = v_["2Y"]*v_["T"*string(Int64(v_["J-1"]))]
                v_["T"*string(J)] = v_["T"*string(J)]-v_["T"*string(Int64(v_["J-2"]))]
            end
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["A0"]
            pbm.A[ig,iv] += Float64(v_["T0"])
            iv = ix_["A2"]
            pbm.A[ig,iv] += Float64(v_["T2"])
            iv = ix_["A4"]
            pbm.A[ig,iv] += Float64(v_["T4"])
            iv = ix_["A6"]
            pbm.A[ig,iv] += Float64(v_["T6"])
            iv = ix_["A8"]
            pbm.A[ig,iv] += Float64(v_["T8"])
            iv = ix_["A10"]
            pbm.A[ig,iv] += Float64(v_["T10"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
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
        pb.x0[ix_["A0"]] = Float64(1.9264e+01)
        pb.x0[ix_["A2"]] = Float64(-1.7302)
        pb.x0[ix_["A4"]] = Float64(4.0794e+01)
        pb.x0[ix_["A6"]] = Float64(8.3021e-01)
        pb.x0[ix_["A8"]] = Float64(3.7090)
        pb.x0[ix_["A10"]] = Float64(-1.7723e-01)
        pb.x0[ix_["K"]] = Float64(10.0)
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
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["12"]):Int64(v_["M"])
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
# LO SOLTN              1.48003482D-04
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-SBR2-RN-8-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
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

