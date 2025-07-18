function DIXMAANP(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DIXMAANP
#    *********
#    A variant on the Dixon-Maany test problem (version I)
# 
#    Source:
#    L. Luksan, C. Matonoha and J. Vlcek  
#    Modified CUTE problems for sparse unconstraoined optimization
#    Technical Report 1081
#    Institute of Computer Science
#    Academy of Science of the Czech Republic
# 
#    (problem 22) based on
# 
#    L.C.W. Dixon and Z. Maany,
#    "A family of test problems with sparse Hessians for unconstrained
#    optimization",
#    TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January 1995.
#               this version Nick Gould, June, 2013
# 
#    classification = "C-COUR2-AN-V-0"
# 
#    M is equal to the third of the number of variables
# 
#       Alternative values for the SIF file parameters:
# IE M                   5              $-PARAMETER n = 15  original value 
# IE M                   30             $-PARAMETER n = 90
# IE M                   100            $-PARAMETER n = 300
# IE M                   500            $-PARAMETER n = 1500
# IE M                   1000           $-PARAMETER n = 3000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DIXMAANP"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling DIXMAANP.")
    end

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
            v_["M"] = Int64(5);  #  SIF file default value
        else
            v_["M"] = Int64(args[1]);
        end
# IE M                   3000           $-PARAMETER n = 9000
        v_["N"] = 3*v_["M"]
        v_["ALPHA"] = 1.0
        v_["BETA"] = 0.26
        v_["GAMMA"] = 0.26
        v_["DELTA"] = 0.26
        v_["K1"] = 2
        v_["K2"] = 1
        v_["K3"] = 1
        v_["K4"] = 2
        v_["RN"] = Float64(v_["N"])
        v_["N-1"] = -1+v_["N"]
        v_["2M"] = v_["M"]+v_["M"]
        v_["1"] = 1
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("GA",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("GB",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("GC",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("GD",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["GA"]] = Float64(-1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(2.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "eSQB", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eSQC", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQB")
            arrset(ielftype,ie,iet_["eSQB"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["2M"])
            v_["I+M"] = I+v_["M"]
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQC")
            arrset(ielftype,ie,iet_["eSQC"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+M"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["I+2M"] = I+v_["2M"]
            ename = "D"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I+2M"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(2.0))
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["I/N"] = v_["RI"]/v_["RN"]
            v_["TMP"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["K1"])
                v_["TMP"] = v_["TMP"]*v_["I/N"]
            end
            v_["AI"] = v_["TMP"]*v_["ALPHA"]
            ig = ig_["GA"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["AI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RI"] = Float64(I)
            v_["I/N"] = v_["RI"]/v_["RN"]
            v_["TMP"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["K2"])
                v_["TMP"] = v_["TMP"]*v_["I/N"]
            end
            v_["BI"] = v_["TMP"]*v_["BETA"]
            ig = ig_["GB"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["BI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["2M"])
            v_["RI"] = Float64(I)
            v_["I/N"] = v_["RI"]/v_["RN"]
            v_["TMP"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["K3"])
                v_["TMP"] = v_["TMP"]*v_["I/N"]
            end
            v_["CI"] = v_["TMP"]*v_["GAMMA"]
            ig = ig_["GC"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["CI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["RI"] = Float64(I)
            v_["I/N"] = v_["RI"]/v_["RN"]
            v_["TMP"] = 1.0
            for J = Int64(v_["1"]):Int64(v_["K4"])
                v_["TMP"] = v_["TMP"]*v_["I/N"]
            end
            v_["DI"] = v_["TMP"]*v_["DELTA"]
            ig = ig_["GD"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(v_["DI"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               1.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-COUR2-AN-V-0"
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

    elseif action == "eSQB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        F1 = EV_[1]*EV_[1]
        F2 = EV_[2]+EV_[2]*EV_[2]
        DF2DY = 1.0+2.0*EV_[2]
        f_   = F1*F2*F2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]*F2*F2
            g_[2] = 2.0*F1*F2*DF2DY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*F2*F2
                H_[1,2] = 4.0*EV_[1]*DF2DY*F2
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0*F1*F2+2.0*F1*DF2DY*DF2DY
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        F1 = EV_[1]*EV_[1]
        F2 = EV_[2]^4
        f_   = F1*F2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]*F2
            g_[2] = 4.0*F1*EV_[2]^3
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*F2
                H_[1,2] = 8.0*EV_[1]*EV_[2]^3
                H_[2,1] = H_[1,2]
                H_[2,2] = 12.0*F1*EV_[2]^2
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

