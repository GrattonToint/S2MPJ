function TWOBARS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Structureal analysis of the simplest two bar scheme.  The structure has
#    the following simple symmetric shape
# 
#                                 *
#                                / \
#                               /   \
#                              /     \
#                            """     """
# 
#    and a force is applied at the top node.  The unknown are the distance
#    of the left and right feet wrt to the projection of the top node and the
#    weight of the bars.
# 
#    Source:
#    an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.
# 
#    SIF input: Ph. Toint, November 1994
# 
#    classification = "C-OOR2-MN-2-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TWOBARS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
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
        ig,ig_,_ = s2mpj_ii("CONS1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CONS1")
        ig,ig_,_ = s2mpj_ii("CONS2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"CONS2")
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
        pbm.gconst[ig_["CONS1"]] = Float64(1.0)
        pbm.gconst[ig_["CONS2"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 0.2
        pb.xupper[ix_["X1"]] = 4.0
        pb.xlower[ix_["X2"]] = 0.1
        pb.xupper[ix_["X2"]] = 1.6
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eOE", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        it,iet_,_ = s2mpj_ii( "eCE1", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        it,iet_,_ = s2mpj_ii( "eCE2", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OBEL"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eOE")
        arrset(ielftype,ie,iet_["eOE"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "COEL1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCE1")
        arrset(ielftype,ie,iet_["eCE1"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "COEL2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCE2")
        arrset(ielftype,ie,iet_["eCE2"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OBEL"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["CONS1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["COEL1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.124))
        ig = ig_["CONS2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["COEL2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.124))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               1.5086379655
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-MN-2-2"
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

    elseif action == "eOE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = 1.0+EV_[2]*EV_[2]
        RA = sqrt(A)
        f_   = EV_[1]*RA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = RA
            g_[2] = EV_[1]*EV_[2]/RA
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = EV_[2]/RA
                H_[2,1] = H_[1,2]
                H_[2,2] = EV_[1]/(A*RA)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCE1"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = 1.0+EV_[2]*EV_[2]
        RA = sqrt(A)
        B = 8.0/EV_[1]
        DB = -8.0/EV_[1]^2
        D2B = 16.0/EV_[1]^3
        C = 1.0/(EV_[1]*EV_[2])
        DCDX = -1.0/(EV_[1]^2*EV_[2])
        DCDY = -1.0/(EV_[2]^2*EV_[1])
        D2CDXX = 2.0/(EV_[1]^3*EV_[2])
        D2CDXY = 1.0/(EV_[1]*EV_[2])^2
        D2CDYY = 2.0/(EV_[1]*EV_[2]^3)
        BC = B+C
        f_   = RA*BC
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = RA*(DB+DCDX)
            g_[2] = EV_[2]*BC/RA+RA*DCDY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = RA*(D2B+D2CDXX)
                H_[1,2] = RA*D2CDXY+EV_[2]*(DB+DCDX)/RA
                H_[2,1] = H_[1,2]
                H_[2,2] = (BC+2.0*EV_[2]*DCDY-EV_[2]*EV_[2]*BC/A)/RA+RA*D2CDYY
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCE2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = 1.0+EV_[2]*EV_[2]
        RA = sqrt(A)
        B = 8.0/EV_[1]
        DB = -8.0/EV_[1]^2
        D2B = 16.0/EV_[1]^3
        C = 1.0/(EV_[1]*EV_[2])
        DCDX = -1.0/(EV_[1]^2*EV_[2])
        DCDY = -1.0/(EV_[2]^2*EV_[1])
        D2CDXX = 2.0/(EV_[1]^3*EV_[2])
        D2CDXY = 1.0/(EV_[1]*EV_[2])^2
        D2CDYY = 2.0/(EV_[1]*EV_[2]^3)
        BC = B-C
        f_   = RA*BC
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = RA*(DB-DCDX)
            g_[2] = EV_[2]*BC/RA-RA*DCDY
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = RA*(D2B-D2CDXX)
                H_[1,2] = -RA*D2CDXY+EV_[2]*(DB-DCDX)/RA
                H_[2,1] = H_[1,2]
                H_[2,2] = (BC-2.0*EV_[2]*DCDY-EV_[2]*EV_[2]*BC/A)/RA-RA*D2CDYY
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

