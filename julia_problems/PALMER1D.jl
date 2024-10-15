function PALMER1D(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER1D
#    *********
# 
#    A linear least squares problem
#    arising from chemical kinetics.
# 
#    model: H-N=N=N TZVP+MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + A12 X**12
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "C-QUR2-RN-7-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER1D"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 35
        v_["1"] = 1
        v_["X1"] = -1.788963
        v_["X2"] = -1.745329
        v_["X3"] = -1.658063
        v_["X4"] = -1.570796
        v_["X5"] = -1.483530
        v_["X6"] = -1.396263
        v_["X7"] = -1.308997
        v_["X8"] = -1.218612
        v_["X9"] = -1.134464
        v_["X10"] = -1.047198
        v_["X11"] = -0.872665
        v_["X12"] = -0.698132
        v_["X13"] = -0.523599
        v_["X14"] = -0.349066
        v_["X15"] = -0.174533
        v_["X16"] = 0.0000000
        v_["X17"] = 1.788963
        v_["X18"] = 1.745329
        v_["X19"] = 1.658063
        v_["X20"] = 1.570796
        v_["X21"] = 1.483530
        v_["X22"] = 1.396263
        v_["X23"] = 1.308997
        v_["X24"] = 1.218612
        v_["X25"] = 1.134464
        v_["X26"] = 1.047198
        v_["X27"] = 0.872665
        v_["X28"] = 0.698132
        v_["X29"] = 0.523599
        v_["X30"] = 0.349066
        v_["X31"] = 0.174533
        v_["X32"] = -1.8762289
        v_["X33"] = -1.8325957
        v_["X34"] = 1.8762289
        v_["X35"] = 1.8325957
        v_["Y1"] = 78.596218
        v_["Y2"] = 65.77963
        v_["Y3"] = 43.96947
        v_["Y4"] = 27.038816
        v_["Y5"] = 14.6126
        v_["Y6"] = 6.2614
        v_["Y7"] = 1.538330
        v_["Y8"] = 0.000000
        v_["Y9"] = 1.188045
        v_["Y10"] = 4.6841
        v_["Y11"] = 16.9321
        v_["Y12"] = 33.6988
        v_["Y13"] = 52.3664
        v_["Y14"] = 70.1630
        v_["Y15"] = 83.4221
        v_["Y16"] = 88.3995
        v_["Y17"] = 78.596218
        v_["Y18"] = 65.77963
        v_["Y19"] = 43.96947
        v_["Y20"] = 27.038816
        v_["Y21"] = 14.6126
        v_["Y22"] = 6.2614
        v_["Y23"] = 1.538330
        v_["Y24"] = 0.000000
        v_["Y25"] = 1.188045
        v_["Y26"] = 4.6841
        v_["Y27"] = 16.9321
        v_["Y28"] = 33.6988
        v_["Y29"] = 52.3664
        v_["Y30"] = 70.1630
        v_["Y31"] = 83.4221
        v_["Y32"] = 108.18086
        v_["Y33"] = 92.733676
        v_["Y34"] = 108.18086
        v_["Y35"] = 92.733676
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
        iv,ix_,_ = s2mpj_ii("A12",ix_)
        arrset(pb.xnames,iv,"A12")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            v_["XQUART"] = v_["XSQR"]*v_["XSQR"]
            v_["X**6"] = v_["XSQR"]*v_["XQUART"]
            v_["X**8"] = v_["XSQR"]*v_["X**6"]
            v_["X**10"] = v_["XSQR"]*v_["X**8"]
            v_["X**12"] = v_["XSQR"]*v_["X**10"]
            v_["X**14"] = v_["XSQR"]*v_["X**12"]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
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
            iv = ix_["A12"]
            pbm.A[ig,iv] += Float64(v_["X**12"])
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
        pb.xlower[ix_["A12"]] = -Inf
        pb.xupper[ix_["A12"]] = +Inf
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
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
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN              0.652673985
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-QUR2-RN-7-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

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

