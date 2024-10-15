function SNAIL(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SNAIL
#    *********
# 
#    A 2D problem featuring a spiraling valley.
#    Dedicated to the city of Namur, whose emblem is a snail.
# 
#    Source:
#    J. Engels, private communication.
# 
#    SIF input: Ph. Toint, May 1990.
# 
#    classification = "C-OUR2-AN-2-0"
# 
#    Problem parameters (CUP > CLOW > 0)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SNAIL"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["CLOW"] = 1.0
        v_["CUP"] = 2.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(10.0)
        pb.x0[ix_["X2"]] = Float64(10.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSPIRAL", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"CL")
        loaset(elftp,it,2,"CU")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSPIRAL")
        arrset(ielftype,ie,iet_["eSPIRAL"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="CL",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["CLOW"]))
        posep = findfirst(x->x=="CU",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["CUP"]))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E"])
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "C-OUR2-AN-2-0"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        return pb, pbm

# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eSPIRAL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = 0.5*(pbm.elpar[iel_][2]+pbm.elpar[iel_][1])
        B = 0.5*(pbm.elpar[iel_][2]-pbm.elpar[iel_][1])
        X2 = EV_[1]*EV_[1]
        Y2 = EV_[2]*EV_[2]
        R2 = X2+Y2
        D = 1.0+R2
        D2 = D*D
        D3 = D2*D
        U = R2/D
        DUDX = (EV_[1]+EV_[1])/D2
        DUDY = (EV_[2]+EV_[2])/D2
        D2UDX2 = 2.0*(D-4.0*X2)/D3
        D2UDY2 = 2.0*(D-4.0*Y2)/D3
        D2UDXY = -8.0*EV_[1]*EV_[2]/D3
        THETA = atan(EV_[2],EV_[1])
        DTDX = -EV_[2]/R2
        DTDY = EV_[1]/R2
        R4 = R2*R2
        D2TDX2 = 2.0*EV_[1]*EV_[2]/R4
        D2TDY2 = -2.0*EV_[2]*EV_[1]/R4
        D2TDXY = (Y2-X2)/R4
        R = sqrt(R2)
        R3 = R*R2
        DRDX = EV_[1]/R
        DRDY = EV_[2]/R
        D2RDX2 = Y2/R3
        D2RDY2 = X2/R3
        D2RDXY = -EV_[1]*EV_[2]/R3
        ARG = R-THETA
        S = B*sin(ARG)
        C = B*cos(ARG)
        DCDX = -S*(DRDX-DTDX)
        DCDY = -S*(DRDY-DTDY)
        D2CDX2 = -C*(DRDX-DTDX)^2-S*(D2RDX2-D2TDX2)
        D2CDY2 = -C*(DRDY-DTDY)^2-S*(D2RDY2-D2TDY2)
        D2CDXY = -C*(DRDX-DTDX)*(DRDY-DTDY)-S*(D2RDXY-D2TDXY)
        V = 1.0+A*R-R*C
        DVDX = A*DRDX-DRDX*C-R*DCDX
        DVDY = A*DRDY-DRDY*C-R*DCDY
        D2VDX2 = A*D2RDX2-D2RDX2*C-2.0*DRDX*DCDX-R*D2CDX2
        D2VDY2 = A*D2RDY2-D2RDY2*C-2.0*DRDY*DCDY-R*D2CDY2
        D2VDXY = A*D2RDXY-D2RDXY*C-DRDX*DCDY-DRDY*DCDX-R*D2CDXY
        f_   = U*V
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DUDX*V+U*DVDX
            g_[2] = DUDY*V+U*DVDY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = D2UDX2*V+2.0*DUDX*DVDX+U*D2VDX2
                H_[1,2] = D2UDXY*V+DUDX*DVDY+DUDY*DVDX+U*D2VDXY
                H_[2,1] = H_[1,2]
                H_[2,2] = D2UDY2*V+2.0*DUDY*DVDY+U*D2VDY2
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

