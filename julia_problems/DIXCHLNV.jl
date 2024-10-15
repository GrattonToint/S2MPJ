function DIXCHLNV(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DIXCHLNV
#    *********
# 
#    A variant of a constrained problem set as a challenge for SQP methods
#    by L.C.W. Dixon at the APMOD91 Conference.
#    The variant from the original problem is that the variables have
#    been constrained to be above 1.0D-15, which in turn allows the
#    constraints to be expressed in terms of logarithms.
# 
#    Source: a modification (by Ph. Toint) of
#    L.C.W. Dixon, personnal communication, Jan 1991.
# 
#    SIF input: Ph. Toint, Feb 1991.
# 
#    classification = "C-SOR2-AN-V-V"
# 
#    Number of variables
#    (variable, but must be even and at least equal to 4)
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER     original value
# IE N                   1000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DIXCHLNV"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   5000           $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["N-3"] = -3+v_["N"]
        v_["N-1"] = -1+v_["N"]
        v_["90.0"] = 90.0
        v_["10.1"] = 10.1
        v_["19.8"] = 19.8
        v_["1/90.0"] = 1.0/v_["90.0"]
        v_["1/10.1"] = 1.0/v_["10.1"]
        v_["1/19.8"] = 1.0/v_["19.8"]
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
        for I = Int64(v_["1"]):Int64(v_["N-3"])
            v_["I+1"] = 1+I
            v_["I+2"] = 2+I
            v_["I+3"] = 3+I
            ig,ig_,_ = s2mpj_ii("A"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(0.01))
            ig,ig_,_ = s2mpj_ii("B"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+3"]))]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(v_["1/90.0"]))
            ig,ig_,_ = s2mpj_ii("D"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+2"]))]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("E"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(v_["1/10.1"]))
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(Int64(v_["I+3"]))]
            pbm.A[ig,iv] += Float64(1.0)
            arrset(pbm.gscale,ig,Float64(v_["1/10.1"]))
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["1/19.8"]))
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("P"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"P"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["N-3"])
            pbm.gconst[ig_["B"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["D"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["E"*string(I)]] = Float64(1.0)
            pbm.gconst[ig_["F"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.0e-8,pb.n)
        pb.xupper = fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        v_["X0"] = 2.0
        for I = Int64(v_["1"]):Int64(v_["2"]):Int64(v_["N-1"])
            pb.x0[ix_["X"*string(I)]] = Float64(v_["X0"])
            v_["1/X0"] = 1.0/v_["X0"]
            v_["I+1"] = 1+I
            pb.x0[ix_["X"*string(Int64(v_["I+1"]))]] = Float64(v_["1/X0"])
            v_["X0"] = 1.0+v_["X0"]
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eS2PR", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        it,iet_,_ = s2mpj_ii( "eLOGV", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "XSQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["N-3"])
            v_["I+1"] = 1+I
            v_["I+3"] = 3+I
            ename = "PR"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eS2PR")
            arrset(ielftype,ie,iet_["eS2PR"])
            vname = "X"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+3"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "LOGX"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eLOGV")
            arrset(ielftype,ie,iet_["eLOGV"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-8),nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N-3"])
            v_["I+2"] = 2+I
            ig = ig_["A"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["B"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["C"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["XSQ"*string(Int64(v_["I+2"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            ig = ig_["D"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["E"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["F"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["PR"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(I)
                ig = ig_["P"*string(I)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["LOGX"*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
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
        pb.pbclass = "C-SOR2-AN-V-V"
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

    elseif action == "eLOGV"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = log(EV_[1])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0/EV_[1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -1.0/(EV_[1]*EV_[1])
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eS2PR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]-1.0)*(EV_[2]-1.0)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]-1.0
            g_[2] = EV_[1]-1.0
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

