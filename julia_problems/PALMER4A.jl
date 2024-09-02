function PALMER4A(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER4A
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#    model: H-N=C=Se TZVP + MP2
#    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6
#                 + B / ( C + X**2 ), B, C nonnegative.
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
# 
#    classification = "SBR2-RN-6-0"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER4A"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

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
        v_["X7"] = -0.741119
        v_["X8"] = -0.698132
        v_["X9"] = -0.523599
        v_["X10"] = -0.349066
        v_["X11"] = -0.174533
        v_["X12"] = 0.0
        v_["X13"] = 0.174533
        v_["X14"] = 0.349066
        v_["X15"] = 0.523599
        v_["X16"] = 0.698132
        v_["X17"] = 0.741119
        v_["X18"] = 0.872665
        v_["X19"] = 1.047198
        v_["X20"] = 1.221730
        v_["X21"] = 1.396263
        v_["X22"] = 1.570796
        v_["X23"] = 1.658063
        v_["Y1"] = 67.27625
        v_["Y2"] = 52.8537
        v_["Y3"] = 30.2718
        v_["Y4"] = 14.9888
        v_["Y5"] = 5.5675
        v_["Y6"] = 0.92603
        v_["Y7"] = 0.0
        v_["Y8"] = 0.085108
        v_["Y9"] = 1.867422
        v_["Y10"] = 5.014768
        v_["Y11"] = 8.263520
        v_["Y12"] = 9.8046208
        v_["Y13"] = 8.263520
        v_["Y14"] = 5.014768
        v_["Y15"] = 1.867422
        v_["Y16"] = 0.085108
        v_["Y17"] = 0.0
        v_["Y18"] = 0.92603
        v_["Y19"] = 5.5675
        v_["Y20"] = 14.9888
        v_["Y21"] = 30.2718
        v_["Y22"] = 52.8537
        v_["Y23"] = 67.27625
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
        iv,ix_,_ = s2mpj_ii("B",ix_)
        arrset(pb.xnames,iv,"B")
        iv,ix_,_ = s2mpj_ii("C",ix_)
        arrset(pb.xnames,iv,"C")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            v_["XQUART"] = v_["XSQR"]*v_["XSQR"]
            v_["XSEXT"] = v_["XQUART"]*v_["XSQR"]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["A0"]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["A2"]
            pbm.A[ig,iv] += Float64(v_["XSQR"])
            iv = ix_["A4"]
            pbm.A[ig,iv] += Float64(v_["XQUART"])
            iv = ix_["A6"]
            pbm.A[ig,iv] += Float64(v_["XSEXT"])
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
        pb.xlower[ix_["B"]] = 0.00001
        pb.xlower[ix_["C"]] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eQUOT", iet_)
        loaset(elftv,it,1,"B")
        loaset(elftv,it,2,"C")
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
            arrset(ielftype, ie, iet_["eQUOT"])
            vname = "B"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="B",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,1.0)
            posev = findfirst(x->x=="C",elftv[ielftype[ie]])
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
# LO SOLTN               4.0606141D-02
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "SBR2-RN-6-0"
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
        DENOM = 1.0/(EV_[2]+pbm.elpar[iel_][1])
        f_   = EV_[1]*DENOM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DENOM
            g_[2] = -EV_[1]*DENOM*DENOM
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -DENOM*DENOM
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[1]*DENOM^3
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

