function ORTHREGF(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ORTHREGF
#    *********
# 
#    An orthogonal regression problem
# 
#    The problem is to fit (orthogonally) an torus to a
#    set of points in 3D space. This set of points is generated by
#    perturbing a first set lying exactly on a predefined torus
#    centered at the origin.
# 
#    Source:
#    M. Gulliksson,
#    "Algorithms for nonlinear Least-squares with Applications to
#    Orthogonal Regression",
#    UMINF-178.90, University of Umea, Sweden, 1990.
# 
#    SIF input: Ph. Toint, June 1990.
#               minor correction by Ph. Shott, Jan 1995.
# 
#    classification = "C-CQOR2-AY-V-V"
# 
#    square root of the number of data points
#    (number of variables = 3 * NPTS**2 + 5 )
# 
#       Alternative values for the SIF file parameters:
# IE NPTS                5              $-PARAMETER n = 80    original value
# IE NPTS                7              $-PARAMETER n = 152
# IE NPTS                10             $-PARAMETER n = 305
# IE NPTS                15             $-PARAMETER n = 680
# IE NPTS                20             $-PARAMETER n = 1205
# IE NPTS                40             $-PARAMETER n = 4805
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 17 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ORTHREGF"

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
            v_["NPTS"] = Int64(5);  #  SIF file default value
        else
            v_["NPTS"] = Int64(args[1]);
        end
        v_["TP4"] = 1.7
        v_["TP5"] = 0.8
        v_["PSEED"] = 237.1531
        v_["PSIZE"] = 0.2
        v_["1"] = 1
        v_["5"] = 5
        v_["PI"] = 3.1415926535
        v_["2PI"] = 2.0*v_["PI"]
        v_["RNPTS"] = Float64(v_["NPTS"])
        v_["ICR0"] = 1.0/v_["RNPTS"]
        v_["INCR"] = v_["ICR0"]*v_["2PI"]
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            v_["I-1"] = -1+I
            v_["RI-1"] = Float64(v_["I-1"])
            v_["THETA1"] = v_["RI-1"]*v_["INCR"]
            v_["ST1"] = sin(v_["THETA1"])
            v_["CT1"] = cos(v_["THETA1"])
            v_["P5CT1"] = v_["TP5"]*v_["CT1"]
            v_["P4P5CT1"] = v_["TP4"]+v_["P5CT1"]
            v_["R3"] = v_["TP5"]*v_["ST1"]
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                v_["J-1"] = -1+J
                v_["RJ-1"] = Float64(v_["J-1"])
                v_["THETA2"] = v_["RJ-1"]*v_["INCR"]
                v_["ST2"] = sin(v_["THETA2"])
                v_["CT2"] = cos(v_["THETA2"])
                v_["R1"] = v_["P4P5CT1"]*v_["CT2"]
                v_["R2"] = v_["P4P5CT1"]*v_["ST2"]
                v_["XSEED"] = v_["THETA2"]*v_["PSEED"]
                v_["SSEED"] = cos(v_["XSEED"])
                v_["PER-1"] = v_["PSIZE"]*v_["SSEED"]
                v_["PERT"] = 1.0+v_["PER-1"]
                v_["XD"*string(I)*","*string(J)] = v_["R1"]*v_["PERT"]
                v_["YD"*string(I)*","*string(J)] = v_["R2"]*v_["PERT"]
                v_["ZD"*string(I)*","*string(J)] = v_["R3"]*v_["PERT"]
            end
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["5"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                iv,ix_,_ = s2mpj_ii("X"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"X"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("Y"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Y"*string(I)*","*string(J))
                iv,ix_,_ = s2mpj_ii("Z"*string(I)*","*string(J),ix_)
                arrset(pb.xnames,iv,"Z"*string(I)*","*string(J))
            end
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                ig,ig_,_ = s2mpj_ii("OX"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("OY"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Y"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("OZ"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Z"*string(I)*","*string(J)]
                pbm.A[ig,iv] += Float64(1.0)
                ig,ig_,_ = s2mpj_ii("A"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"A"*string(I)*","*string(J))
            end
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
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                pbm.gconst[ig_["OX"*string(I)*","*string(J)]]  = (
                      Float64(v_["XD"*string(I)*","*string(J)]))
                pbm.gconst[ig_["OY"*string(I)*","*string(J)]]  = (
                      Float64(v_["YD"*string(I)*","*string(J)]))
                pbm.gconst[ig_["OZ"*string(I)*","*string(J)]]  = (
                      Float64(v_["ZD"*string(I)*","*string(J)]))
            end
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["P4"]] = 0.001
        pb.xlower[ix_["P5"]] = 0.001
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"P1")
            pb.x0[ix_["P1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["P1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"P2")
            pb.x0[ix_["P2"]] = Float64(0.0)
        else
            pb.y0[findfirst(x->x==ig_["P2"],pbm.congrps)] = Float64(0.0)
        end
        if haskey(ix_,"P3")
            pb.x0[ix_["P3"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["P3"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"P4")
            pb.x0[ix_["P4"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["P4"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"P5")
            pb.x0[ix_["P5"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["P5"],pbm.congrps)] = Float64(0.5)
        end
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                if haskey(ix_,"X"*string(I)*","*string(J))
                    pb.x0[ix_["X"*string(I)*","*string(J)]]  = (
                          Float64(v_["XD"*string(I)*","*string(J)]))
                else
                    pb.y0[findfirst(x->x==ig_["X"*string(I)*","*string(J)],pbm.congrps)]  = (
                          Float64(v_["XD"*string(I)*","*string(J)]))
                end
                if haskey(ix_,"Y"*string(I)*","*string(J))
                    pb.x0[ix_["Y"*string(I)*","*string(J)]]  = (
                          Float64(v_["YD"*string(I)*","*string(J)]))
                else
                    pb.y0[findfirst(x->x==ig_["Y"*string(I)*","*string(J)],pbm.congrps)]  = (
                          Float64(v_["YD"*string(I)*","*string(J)]))
                end
                if haskey(ix_,"Z"*string(I)*","*string(J))
                    pb.x0[ix_["Z"*string(I)*","*string(J)]]  = (
                          Float64(v_["ZD"*string(I)*","*string(J)]))
                else
                    pb.y0[findfirst(x->x==ig_["Z"*string(I)*","*string(J)],pbm.congrps)]  = (
                          Float64(v_["ZD"*string(I)*","*string(J)]))
                end
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eTA", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        loaset(elftv,it,3,"A")
        loaset(elftv,it,4,"B")
        loaset(elftv,it,5,"C")
        it,iet_,_ = s2mpj_ii( "eISQ", iet_)
        loaset(elftv,it,1,"Z")
        loaset(elftv,it,2,"P")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"XX")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["NPTS"])
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                ename = "EA"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eTA")
                arrset(ielftype,ie,iet_["eTA"])
                vname = "X"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "Y"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "P1"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="A",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "P2"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="B",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "P4"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="C",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "EB"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eISQ")
                arrset(ielftype,ie,iet_["eISQ"])
                vname = "Z"*string(I)*","*string(J)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                vname = "P3"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="P",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
                ename = "EC"*string(I)*","*string(J)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "P5"
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
                posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
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
            for J = Int64(v_["1"]):Int64(v_["NPTS"])
                ig = ig_["OX"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["OY"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["OZ"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gL2")
                ig = ig_["A"*string(I)*","*string(J)]
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EA"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,1.)
                posel = posel+1
                loaset(pbm.grelt,ig,posel,ie_["EB"*string(I)*","*string(J)])
                loaset(pbm.grelw,ig,posel, 1.)
                posel = length(pbm.grelt[ig])+1
                loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)*","*string(J)])
                arrset(nlc,length(nlc)+1,ig)
                loaset(pbm.grelw,ig,posel,Float64(-1.0))
            end
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN(5)            0.990089426
# LO SOLTN(7)            1.315031322
# LO SOLTN(10)           4.515848902
# LO SOLTN(15)           9.185538338
# LO SOLTN(20)           16.20054380
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
        pb.pbclass = "C-CQOR2-AY-V-V"
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
        CCSQ = IV_[3]*IV_[3]
        CCCB = CCSQ*IV_[3]
        XXYY = IV_[1]*IV_[1]+IV_[2]*IV_[2]
        T = XXYY/CCSQ
        DTDX = 2.0*IV_[1]/CCSQ
        DTDY = 2.0*IV_[2]/CCSQ
        DTDC = -2.0*XXYY/CCCB
        D2TDX2 = 2.0/CCSQ
        D2TDY2 = 2.0/CCSQ
        D2TDC2 = 6.0*XXYY/(CCSQ*CCSQ)
        D2TDXC = -4.0*IV_[1]/CCCB
        D2TDYC = -4.0*IV_[2]/CCCB
        S = sqrt(T)
        R = 0.5/S
        DSDX = R*DTDX
        DSDY = R*DTDY
        DSDC = R*DTDC
        D2SDX2 = R*(D2TDX2-0.5*DTDX*DTDX/T)
        D2SDY2 = R*(D2TDY2-0.5*DTDY*DTDY/T)
        D2SDC2 = R*(D2TDC2-0.5*DTDC*DTDC/T)
        D2SDXY = -0.5*DTDX*DSDY/T
        D2SDXC = R*(D2TDXC-0.5*DTDX*DTDC/T)
        D2SDYC = R*(D2TDYC-0.5*DTDY*DTDC/T)
        SS = S-1.0
        SPS = SS+SS
        f_   = SS*SS
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = SPS*DSDX
            g_[2] = SPS*DSDY
            g_[3] = SPS*DSDC
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = SPS*D2SDX2+2.0*DSDX*DSDX
                H_[1,2] = SPS*D2SDXY+2.0*DSDX*DSDY
                H_[2,1] = H_[1,2]
                H_[1,3] = SPS*D2SDXC+2.0*DSDX*DSDC
                H_[3,1] = H_[1,3]
                H_[2,2] = SPS*D2SDY2+2.0*DSDY*DSDY
                H_[2,3] = SPS*D2SDYC+2.0*DSDY*DSDC
                H_[3,2] = H_[2,3]
                H_[3,3] = SPS*D2SDC2+2.0*DSDC*DSDC
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

    elseif action == "eISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

