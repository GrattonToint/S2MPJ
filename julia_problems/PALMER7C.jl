function PALMER7C(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER7C
#    *********
# 
#    A linear least squares problem
#    arising from chemical kinetics.
# 
#     model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
#                 A10 X**10 + A12 X**12 + A14 X**14
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1992.
# 
#    classification = "QUR2-RN-8-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER7C"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "PALMER7C"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 24
        v_["1"] = 1
        v_["12"] = 12
        v_["X12"] = 0.000000
        v_["X13"] = 0.139626
        v_["X14"] = 0.261799
        v_["X15"] = 0.436332
        v_["X16"] = 0.565245
        v_["X17"] = 0.512942
        v_["X18"] = 0.610865
        v_["X19"] = 0.785398
        v_["X20"] = 0.959931
        v_["X21"] = 1.134464
        v_["X22"] = 1.308997
        v_["X23"] = 1.483530
        v_["X24"] = 1.658063
        v_["Y12"] = 4.419446
        v_["Y13"] = 3.564931
        v_["Y14"] = 2.139067
        v_["Y15"] = 0.404686
        v_["Y16"] = 0.000000
        v_["Y17"] = 0.035152
        v_["Y18"] = 0.146813
        v_["Y19"] = 2.718058
        v_["Y20"] = 9.474417
        v_["Y21"] = 26.132221
        v_["Y22"] = 41.451561
        v_["Y23"] = 72.283164
        v_["Y24"] = 117.630959
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
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
        iv,ix_,_ = s2mpj_ii("A14",ix_)
        arrset(pb.xnames,iv,"A14")
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
            iv = ix_["A14"]
            pbm.A[ig,iv] += Float64(v_["X**14"])
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
        pb.xlower[ix_["A12"]] = -Inf
        pb.xupper[ix_["A12"]] = +Inf
        pb.xlower[ix_["A14"]] = -Inf
        pb.xupper[ix_["A14"]] = +Inf
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
        for I = Int64(v_["12"]):Int64(v_["M"])
            ig = ig_["O"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN              5.0310687D-02
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "QUR2-RN-8-0"
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

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2mpj_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

