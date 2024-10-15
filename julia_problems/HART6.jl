function HART6(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HART6
#    *********
# 
#    Source: Hartman problem 6 in
#    L. C. W. Dixon and G. P. Szego (Eds.)
#    Towards Global Optimization
#    North Holland, 1975.
#    Paper 9, page 163.
# 
#    SIF input: A.R. Conn May 1995
# 
#    classification = "C-OBR2-AN-6-0"
# 
#    Number of variables - constraints
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HART6"

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
        v_["ONE"] = 1
        v_["NN"] = 6
        v_["L"] = 4
        v_["C1"] = 1.0
        v_["C2"] = 1.2
        v_["C3"] = 3.0
        v_["C4"] = 3.2
        v_["A1,1"] = 10.0
        v_["A2,1"] = 0.05
        v_["A3,1"] = 3.0
        v_["A4,1"] = 17.0
        v_["A1,2"] = 0.05
        v_["A2,2"] = 10.0
        v_["A3,2"] = 3.5
        v_["A4,2"] = 8.0
        v_["A1,3"] = 17.0
        v_["A2,3"] = 17.0
        v_["A3,3"] = 1.7
        v_["A4,3"] = 0.05
        v_["A1,4"] = 3.5
        v_["A2,4"] = 0.1
        v_["A3,4"] = 10.0
        v_["A4,4"] = 10.0
        v_["A1,5"] = 1.7
        v_["A2,5"] = 8.0
        v_["A3,5"] = 17.0
        v_["A4,5"] = 0.1
        v_["A1,6"] = 8.0
        v_["A2,6"] = 14.0
        v_["A3,6"] = 8.0
        v_["A4,6"] = 14.0
        v_["P1,1"] = 0.1312
        v_["P2,1"] = 0.2329
        v_["P3,1"] = 0.2348
        v_["P4,1"] = 0.4047
        v_["P1,2"] = 0.1696
        v_["P2,2"] = 0.4135
        v_["P3,2"] = 0.1451
        v_["P4,2"] = 0.8828
        v_["P1,3"] = 0.5569
        v_["P2,3"] = 0.8307
        v_["P3,3"] = 0.3522
        v_["P4,3"] = 0.8732
        v_["P1,4"] = 0.0124
        v_["P2,4"] = 0.3736
        v_["P3,4"] = 0.2883
        v_["P4,4"] = 0.5743
        v_["P1,5"] = 0.8283
        v_["P2,5"] = 0.1004
        v_["P3,5"] = 0.3047
        v_["P4,5"] = 0.1091
        v_["P1,6"] = 0.5886
        v_["P2,6"] = 0.9991
        v_["P3,6"] = 0.6650
        v_["P4,6"] = 0.0381
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
        for I = Int64(v_["1"]):Int64(v_["L"])
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(-1.0))
        end
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.0,pb.n)
        pb.xupper = fill(1.0,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.2),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V1")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PIJ")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["L"])
            for J = Int64(v_["1"]):Int64(v_["NN"])
                ename = "E"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "X"*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),Float64(1.0),Float64(0.2))
                posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                posep = findfirst(x->x=="PIJ",elftp[ielftype[ie]])
                loaset(pbm.elpar,ie,posep,Float64(v_["P"*string(I)*","*string(J)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gNEXP",igt_)
        it,igt_,_ = s2mpj_ii("gNEXP",igt_)
        grftp = Vector{Vector{String}}()
        loaset(grftp,it,1,"CI")
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["L"])
            ig = ig_["OBJ"*string(I)]
            arrset(pbm.grftype,ig,"gNEXP")
            for J = Int64(v_["1"]):Int64(v_["NN"])
                ig = ig_["OBJ"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["E"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel,Float64(v_["A"*string(I)*","*string(J)]))
            end
            ig = ig_["OBJ"*string(I)]
            posgp = findfirst(x->x=="CI",grftp[igt_[pbm.grftype[ig]]])
            loaset(pbm.grpar,ig,posgp,Float64(v_["C"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -3.32288689158
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OBR2-AN-6-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]-pbm.elpar[iel_][1])*(EV_[1]-pbm.elpar[iel_][1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*(EV_[1]-pbm.elpar[iel_][1])
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gNEXP"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= pbm.grpar[igr_][1]*exp(-GVAR_)
        if nargout>1
            g_ = -pbm.grpar[igr_][1]*exp(-GVAR_)
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = pbm.grpar[igr_][1]*exp(-GVAR_)
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

