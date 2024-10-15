function DIXMAANM1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : DIXMAANM1
#    *********
#    A variant on the Dixon-Maany test problem (version I) but removing 
#    elements/groups of type 2 since the parameter beta=0
# 
#    Source:
#    L. Luksan, C. Matonoha and J. Vlcek  
#    Modified CUTE problems for sparse unconstraoined optimization
#    Technical Report 1081
#    Institute of Computer Science
#    Academy of Science of the Czech Republic
# 
#    (problem 19) based on
# 
#    L.C.W. Dixon and Z. Maany,
#    "A family of test problems with sparse Hessians for unconstrained
#    optimization",
#    TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
# 
#    SIF input: Ph. Toint, Dec 1989.
#               correction by Ph. Shott, January 1995.
#               this version Nick Gould, June, 2013
#               update Nick Gould, August 2022, to remove beta=0 terms.
# 
#    classification = "C-OUR2-AN-V-0"
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
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "DIXMAANM1"

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
        v_["BETA"] = 0.0
        v_["GAMMA"] = 0.125
        v_["DELTA"] = 0.125
        v_["K1"] = 2
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("GA",ig_)
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
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(2.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
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
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-V-0"
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

