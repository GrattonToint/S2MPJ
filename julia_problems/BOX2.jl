function BOX2(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : BOX2
#    *********
#    Box problem in 2 variables, obtained by fixing X3 = 1 in BOX3.
# 
#    Source: Problem 11 in
#    A.R. Buckley,
#    "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "SXR2-AN-3-0"
# 
#    This function  is a nonlinear least squares with 10 groups.  Each
#    group has 2 nonlinear elements of exponential type.
# 
#    Number of groups
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "BOX2"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 10
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["MRI"] = -1.0*v_["RI"]
            v_["MTI"] = -0.1*v_["RI"]
            v_["EMTI"] = exp(v_["MTI"])
            v_["EMRI"] = exp(v_["MRI"])
            v_["MEMTI"] = -1.0*v_["EMTI"]
            v_["COEFF"] = v_["MEMTI"]+v_["EMRI"]
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X3"]
            pbm.A[ig,iv] += Float64(v_["COEFF"])
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X3"]] = 1.0
        pb.xupper[ix_["X3"]] = 1.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(0.0)
        pb.x0[ix_["X2"]] = Float64(10.0)
        pb.x0[ix_["X3"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXPT", iet_)
        loaset(elftv,it,1,"V")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"T")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["MTI"] = -0.1*v_["RI"]
            ename = "A"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eEXPT")
                arrset(ielftype,ie,iet_["eEXPT"])
            end
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["MTI"]))
            ename = "B"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eEXPT")
                arrset(ielftype,ie,iet_["eEXPT"])
            end
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="T",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["MTI"]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for ig in 1:ngrp
            arrset(pbm.grftype,ig,"gL2")
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "SXR2-AN-3-0"
        return pb, pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eEXPT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPA = exp(pbm.elpar[iel_][1]*EV_[1])
        f_   = EXPA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*EXPA
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*EXPA
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

