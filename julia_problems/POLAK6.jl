function POLAK6(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POLAK6
#    *********
# 
#    A nonlinear minmax problem in four variables. This is a variation
#    on problem ROSENMMX.
# 
#    Source: 
#    E. Polak, D.H. Mayne and J.E. Higgins,
#    "Superlinearly convergent algorithm for min-max problems"
#    JOTA 69, pp. 407-439, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-LOR2-AN-5-4"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "POLAK6"

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
        iv,ix_,_ = s2mpj_ii("X3",ix_)
        arrset(pb.xnames,iv,"X3")
        iv,ix_,_ = s2mpj_ii("X4",ix_)
        arrset(pb.xnames,iv,"X4")
        iv,ix_,_ = s2mpj_ii("U",ix_)
        arrset(pb.xnames,iv,"U")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("F1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"F1")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-5.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-5.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-21.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(7.0)
        ig,ig_,_ = s2mpj_ii("F2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"F2")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-15.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-11.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig,ig_,_ = s2mpj_ii("F3",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"F3")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-15.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-5.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-21.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-3.0)
        ig,ig_,_ = s2mpj_ii("F4",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"F4")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(15.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-15.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-21.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-3.0)
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
        pbm.gconst[ig_["F2"]] = Float64(80.0)
        pbm.gconst[ig_["F3"]] = Float64(100.0)
        pbm.gconst[ig_["F4"]] = Float64(50.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"XX")
        it,iet_,_ = s2mpj_ii( "eEL42", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        it,iet_,_ = s2mpj_ii( "eEL442", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        loaset(elftv,it,3,"ZZ")
        it,iet_,_ = s2mpj_ii( "eEL4", iet_)
        loaset(elftv,it,1,"XX")
        it,iet_,_ = s2mpj_ii( "eEL44", iet_)
        loaset(elftv,it,1,"XX")
        loaset(elftv,it,2,"YY")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL42")
        arrset(ielftype,ie,iet_["eEL42"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL442")
        arrset(ielftype,ie,iet_["eEL442"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="ZZ",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL4")
        arrset(ielftype,ie,iet_["eEL4"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL44")
        arrset(ielftype,ie,iet_["eEL44"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="YY",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X3SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "X4SQ"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="XX",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["F1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4SQ"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(5.0))
        ig = ig_["F2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(12.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4SQ"])
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-5.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(15.0))
        ig = ig_["F3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel,Float64(21.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(12.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4SQ"])
        loaset(pbm.grelw,ig,posel,Float64(21.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(15.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(5.0))
        ig = ig_["F4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel,Float64(11.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["X3SQ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(12.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["X4SQ"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-15.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel,Float64(15.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution at ( 0, 1, 2, -1 )
# LO SOLTN               -44.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-LOR2-AN-5-4"
        pb.x0          = zeros(Float64,pb.n)
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

    elseif action == "eEL42"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[2]+1.0
        A = EV_[1]-B^4
        f_   = A*A
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*A
            g_[2] = -8.0*A*B^3
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0
                H_[1,2] = -8.0*B^3
                H_[2,1] = H_[1,2]
                H_[2,2] = 32.0*B^6-24.0*A*B^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEL442"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[3]+1.0
        C = EV_[2]-B^4
        DCDZ = -4.0*B^3
        D2CDZZ = -12.0*B^2
        A = EV_[1]-C^4
        DADY = -4.0*C^3
        DADZ = DADY*DCDZ
        D2ADYY = -12.0*C^2
        D2ADYZ = D2ADYY*DCDZ
        D2ADZZ = D2ADYZ*DCDZ+DADY*D2CDZZ
        f_   = A*A
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*A
            g_[2] = 2.0*A*DADY
            g_[3] = 2.0*A*DADZ
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = 2.0
                H_[1,2] = 2.0*DADY
                H_[2,1] = H_[1,2]
                H_[1,3] = 2.0*DADZ
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*(DADY*DADY+A*D2ADYY)
                H_[2,3] = 2.0*(DADZ*DADY+A*D2ADYZ)
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*(DADZ*DADZ+A*D2ADZZ)
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEL4"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[1]+1.0
        f_   = B^4
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 4.0*B^3
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 12.0*B^2
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEL44"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        B = EV_[2]+1.0
        A = EV_[1]-B^4
        f_   = A^4
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 4.0*A^3
            g_[2] = -16.0*(A*B)^3
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 12.0*A^2
                H_[1,2] = -48.0*A^2*B^3
                H_[2,1] = H_[1,2]
                H_[2,2] = -48.0*(A*B)^2*(A-4.0*B^4)
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

