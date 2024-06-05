function HAIRY(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HAIRY
#    *********
# 
#    A hairy problem in two variables.  The surface defined by
#    this function has a large number of relatively sharp hills between
#    which a valley leads to the minimizer.
#    This problem contains a large number of saddle points.
# 
#    Dedicated to Meret Oppenheim, creator of the "furry cup" (1936).
# 
#    Source:
#    Ph. Toint, private communication,
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "OUR2-AY-2-0"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HAIRY"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "HAIRY"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["HLENGTH"] = 30.0
        v_["CSLOPE"] = 100.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("X1",ix_)
        arrset(pb.xnames,iv,"X1")
        iv,ix_,_ = s2mpj_ii("X2",ix_)
        arrset(pb.xnames,iv,"X2")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("FURCUP",ig_)
        arrset(gtype,ig,"<>")
        #%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = length(ix_)
        ngrp   = length(ig_)
        pbm.objgrps = collect(1:ngrp)
        pb.m        = 0
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.x0[ix_["X1"]] = Float64(-5.0)
        pb.x0[ix_["X2"]] = Float64(-7.0)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eFUR", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"DENS")
        it,iet_,_ = s2mpj_ii( "eDCUP", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftp,it,1,"SMOOTH")
        it,iet_,_ = s2mpj_ii( "en1CUP", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftp,it,1,"SMOOTH")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "HAIR"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eFUR")
        arrset(ielftype, ie, iet_["eFUR"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="DENS",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(7.0))
        ename = "DBOWL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eDCUP")
        arrset(ielftype, ie, iet_["eDCUP"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="SMOOTH",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.01))
        ename = "1BOWL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"en1CUP")
        arrset(ielftype, ie, iet_["en1CUP"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="SMOOTH",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.01))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["FURCUP"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["HAIR"])
        loaset(pbm.grelw,ig,posel,Float64(v_["HLENGTH"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["DBOWL"])
        loaset(pbm.grelw,ig,posel,Float64(v_["CSLOPE"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["1BOWL"])
        loaset(pbm.grelw,ig,posel,Float64(v_["CSLOPE"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               20.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = "OUR2-AY-2-0"
        return pb, pbm
# **********************
#  SET UP THE FUNCTION *
#  AND RANGE ROUTINES  *
# **********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eFUR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DV1 = pbm.elpar[iel_][1]*EV_[1]
        DV2 = pbm.elpar[iel_][1]*EV_[2]
        TDV1 = DV1+DV1
        TDV2 = DV2+DV2
        TDL2 = 2.0*pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        S1SQ = sin(DV1)^2
        C2SQ = cos(DV2)^2
        STDV1 = sin(TDV1)
        STDV2 = sin(TDV2)
        f_   = S1SQ*C2SQ
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][1]*STDV1*C2SQ
            g_[2] = -pbm.elpar[iel_][1]*S1SQ*STDV2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = TDL2*cos(TDV1)*C2SQ
                H_[1,2] = -pbm.elpar[iel_][1]*pbm.elpar[iel_][1]*STDV1*STDV2
                H_[2,1] = H_[1,2]
                H_[2,2] = -TDL2*S1SQ*cos(TDV2)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eDCUP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        VSQ = IV_[1]*IV_[1]
        ARG = pbm.elpar[iel_][1]+VSQ
        SQARG = sqrt(ARG)
        DEN = 1.0/SQARG
        f_   = SQARG
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]*DEN
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = (1.0-VSQ/ARG)*DEN
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

    elseif action == "en1CUP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        VSQ = EV_[1]*EV_[1]
        ARG = pbm.elpar[iel_][1]+VSQ
        SQARG = sqrt(ARG)
        DEN = 1.0/SQARG
        f_   = SQARG
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[1]*DEN
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = (1.0-VSQ/ARG)*DEN
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

