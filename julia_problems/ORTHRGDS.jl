function ORTHRGDS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHRGDS
#    *********
# 
#    An orthogonal regression problem.
# 
#    The problem is to fit (orthogonally) a circle to a set of points
#    in the plane. This set of points is generated by perturbing a
#    first set lying exactly on a predefined circle centered at the
#    origin. Each point mapped to the center of the circle has a constraint 
#    gradient of zero, making the Jacobian of the constraints singular.
#    Certain data points (which serve to define the initial estimate
#    of the unknowns) are modified so that they will converge to
#    center of the circle, allowing us to study ill-conditioned
#    Jacobians.
#    With 10 singular points, the final objective is SOLTN = 118.037
# 
#    Source: adapted from:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, Mar 1991 and T. Plantenga, May 1993.
# 
#    classification = "C-QOR2-AY-V-V"
# 
#    Number of data points
#    (number of variables = 2 NPTS + 3 )
# 
#       Alternative values for the SIF file parameters:
# IE NPTS                10             $-PARAMETER n = 23
# IE NPTS                50             $-PARAMETER n = 103
# IE NPTS                76             $-PARAMETER n = 155    original value
# IE NPTS                100            $-PARAMETER n = 203
# IE NPTS                250            $-PARAMETER n = 503
# IE NPTS                500            $-PARAMETER n = 1003
# IE NPTS                2500           $-PARAMETER n = 5003
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ORTHRGDS"

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
            v_["NPTS"] = Int64(20);  #  SIF file default value
        else
            v_["NPTS"] = Int64(args[1]);
        end
# IE NPTS                5000           $-PARAMETER n = 10003
        v_["TZ3"] = 1.7
        v_["PSEED"] = 237.1531
        v_["PSIZE"] = 0.2
        v_["1"] = 1
        v_["0"] = 0
        v_["TDPulo"] = 5
        v_["TDPuhi"] = 14
        v_["PI"] = 3.1415926535
        v_["2PI"] = 2.0*v_["PI"]
        v_["RNPTS"] = Float64(v_["NPTS"])
        v_["ICR0"] = 1.0/v_["RNPTS"]
        v_["INCR"] = v_["ICR0"]*v_["2PI"]
        v_["Z3SQ"] = v_["TZ3"]*v_["TZ3"]
        v_["1+TZ3SQ"] = 1.0+v_["Z3SQ"]
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["THETA"] = v_["RI-1"]*v_["INCR"]
            v_["ST"] = sin(v_["THETA"])
            v_["CT"] = cos(v_["THETA"])
            v_["FACT"] = v_["1+TZ3SQ"]+v_["CT"]
            v_["R1"] = v_["FACT"]*v_["CT"]
            v_["R2"] = v_["FACT"]*v_["ST"]
            v_["XSEED"] = v_["THETA"]*v_["PSEED"]
            v_["SSEED"] = cos(v_["XSEED"])
            v_["PER-1"] = v_["PSIZE"]*v_["SSEED"]
            v_["PERT"] = 1.0+v_["PER-1"]
            v_["XD"*string(I)] = v_["R1"]*v_["PERT"]
            v_["YD"*string(I)] = v_["R2"]*v_["PERT"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("Z1",ix_)
        arrset(pb.xnames,iv,"Z1")
        iv,ix_,_ = s2mpj_ii("Z2",ix_)
        arrset(pb.xnames,iv,"Z2")
        iv,ix_,_ = s2mpj_ii("Z3",ix_)
        arrset(pb.xnames,iv,"Z3")
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ig,ig_,_ = s2mpj_ii("OX"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("OY"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
            ig,ig_,_ = s2mpj_ii("E"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            pbm.gconst[ig_["OX"*string(I)]] = Float64(v_["XD"*string(I)])
            pbm.gconst[ig_["OY"*string(I)]] = Float64(v_["YD"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"Z1")
            pb.x0[ix_["Z1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["Z1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"Z2")
            pb.x0[ix_["Z2"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["Z2"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"Z3")
            pb.x0[ix_["Z3"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["Z3"],pbm.congrps)] = Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(v_["XD"*string(I)])
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)]  = (
                      Float64(v_["XD"*string(I)]))
            end
            if haskey(ix_,"Y"*string(I))
                pb.x0[ix_["Y"*string(I)]] = Float64(v_["YD"*string(I)])
            else
                pb.y0[findfirst(x->x==ig_["Y"*string(I)],pbm.congrps)]  = (
                      Float64(v_["YD"*string(I)]))
            end
        end
        for I = Int64(v_["TDPulo"]):Int64(v_["TDPuhi"])
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(1.8)
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)] = Float64(1.8)
            end
            if haskey(ix_,"Y"*string(I))
                pb.x0[ix_["Y"*string(I)]] = Float64(1.0)
            else
                pb.y0[findfirst(x->x==ig_["Y"*string(I)],pbm.congrps)] = Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eTA", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"ZA")
        loaset(elftv,it,4,"ZB")
        it,iet_,_ = s2mpj_ii( "eTB", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        loaset(elftv,it,3,"ZA")
        loaset(elftv,it,4,"ZB")
        loaset(elftv,it,5,"ZC")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ename = "EA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eTA")
            arrset(ielftype,ie,iet_["eTA"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="ZA",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="ZB",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eTB")
            arrset(ielftype,ie,iet_["eTB"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="ZA",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="ZB",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="ZC",elftv[ielftype[ie]])
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
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            ig = ig_["OX"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["OY"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["E"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(10)           3.412121065
# LO SOLTN(50)           15.59042181
# LO SOLTN(250)          76.10435792
# LO SOLTN(500)          151.2351183
# LO SOLTN(2500)         ???
# LO SOLTN(5000)         ???
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
        pb.pbclass = "C-QOR2-AY-V-V"
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

    elseif action == "eTA"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,4)
        IV_ =  zeros(Float64,2)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        T = IV_[1]*IV_[1]+IV_[2]*IV_[2]
        f_   = T*T
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 4.0*T*IV_[1]
            g_[2] = 4.0*T*IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 4.0*(T+2.0*IV_[1]*IV_[1])
                H_[1,2] = 8.0*IV_[1]*IV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 4.0*(T+2.0*IV_[2]*IV_[2])
                H_ = U_'*H_*U_
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eTB"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,3,5)
        IV_ =  zeros(Float64,3)
        U_[1,1] = U_[1,1]+1
        U_[1,3] = U_[1,3]-1
        U_[2,2] = U_[2,2]+1
        U_[2,4] = U_[2,4]-1
        U_[3,5] = U_[3,5]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        IV_[3] = dot(U_[3,:],EV_)
        T = IV_[1]*IV_[1]+IV_[2]*IV_[2]
        ZZSQ = IV_[3]*IV_[3]
        T1 = 1.0+ZZSQ
        T1SQ = T1*T1
        f_   = T*T1SQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*IV_[1]*T1SQ
            g_[2] = 2.0*IV_[2]*T1SQ
            g_[3] = 4.0*T*T1*IV_[3]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0*T1SQ
                H_[1,3] = 8.0*IV_[1]*T1*IV_[3]
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*T1SQ
                H_[2,3] = 8.0*IV_[2]*T1*IV_[3]
                H_[3,2] = H_[2,3]
                H_[3,3] = 4.0*T*(2.0*ZZSQ+T1)
                H_ = U_'*H_*U_
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

