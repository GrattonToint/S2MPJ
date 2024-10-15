function SYNTHES2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : SYNTHES2
#    *********
# 
#    Source: Test problem 2 (Synthesis of processing system) in 
#    M. Duran & I.E. Grossmann,
#    "An outer approximation algorithm for a class of mixed integer nonlinear
#     programs", Mathematical Programming 36, pp. 307-339, 1986.
# 
#    SIF input: S. Leyffer, October 1997
# 
#    classification = "C-OOR2-AN-11-14"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "SYNTHES2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["1"] = 1
        v_["5"] = 5
        v_["6"] = 6
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["6"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["5"])
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Y1"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["Y2"]
        pbm.A[ig,iv] += Float64(8.0)
        iv = ix_["Y3"]
        pbm.A[ig,iv] += Float64(6.0)
        iv = ix_["Y4"]
        pbm.A[ig,iv] += Float64(10.0)
        iv = ix_["Y5"]
        pbm.A[ig,iv] += Float64(6.0)
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-10.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-15.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-15.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(15.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(5.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-20.0)
        ig,ig_,_ = s2mpj_ii("N1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"N1")
        ig,ig_,_ = s2mpj_ii("N2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"N2")
        iv = ix_["Y1"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("N3",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"N3")
        iv = ix_["Y2"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("L1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L1")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.25)
        iv = ix_["Y3"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("L2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L2")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Y4"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("L3",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L3")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["Y5"]
        pbm.A[ig,iv] += Float64(-10.0)
        ig,ig_,_ = s2mpj_ii("L4",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L4")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-2.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(2.0)
        ig,ig_,_ = s2mpj_ii("L5",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L5")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(-0.75)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(2.0)
        ig,ig_,_ = s2mpj_ii("L6",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L6")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("L7",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L7")
        iv = ix_["X3"]
        pbm.A[ig,iv] += Float64(2.0)
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-1.0)
        iv = ix_["X6"]
        pbm.A[ig,iv] += Float64(-2.0)
        ig,ig_,_ = s2mpj_ii("L8",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L8")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(-0.5)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("L9",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L9")
        iv = ix_["X4"]
        pbm.A[ig,iv] += Float64(0.2)
        iv = ix_["X5"]
        pbm.A[ig,iv] += Float64(-1.0)
        ig,ig_,_ = s2mpj_ii("L10",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"L10")
        iv = ix_["Y1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Y2"]
        pbm.A[ig,iv] += Float64(1.0)
        ig,ig_,_ = s2mpj_ii("L11",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"L11")
        iv = ix_["Y4"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["Y5"]
        pbm.A[ig,iv] += Float64(1.0)
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
        pbm.gconst[ig_["OBJ"]] = Float64(-140.0)
        pbm.gconst[ig_["N2"]] = Float64(1.0)
        pbm.gconst[ig_["N3"]] = Float64(1.0)
        pbm.gconst[ig_["L10"]] = Float64(1.0)
        pbm.gconst[ig_["L11"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        arrset(grange,ig_["L11"],Float64(-1.0))
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xupper[ix_["X1"]] = 2.0
        pb.xupper[ix_["X2"]] = 2.0
        pb.xupper[ix_["X3"]] = 2.0
        pb.xupper[ix_["X6"]] = 3.0
        for I = Int64(v_["1"]):Int64(v_["5"])
            pb.xupper[ix_["Y"*string(I)]] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eLOGSUM", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eEXPA", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"A")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "LOGX4X5"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eLOGSUM")
        arrset(ielftype,ie,iet_["eLOGSUM"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EXPX1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEXPA")
        arrset(ielftype,ie,iet_["eEXPA"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="A",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "EXPX2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEXPA")
        arrset(ielftype,ie,iet_["eEXPA"])
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="A",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.2))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EXPX1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EXPX2"])
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["LOGX4X5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-60.0))
        ig = ig_["N1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["LOGX4X5"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["N2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EXPX1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        ig = ig_["N3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EXPX2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[1:pb.nle] = grange[legrps]
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AN-11-14"
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

    elseif action == "eLOGSUM"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = log(EV_[1]+EV_[2]+1.0)
        DX = 1.0/(EV_[1]+EV_[2]+1.0)
        DXDX = -1.0/(EV_[1]+EV_[2]+1.0)^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DX
            g_[2] = DX
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = DXDX
                H_[1,2] = DXDX
                H_[2,1] = H_[1,2]
                H_[2,2] = DXDX
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eEXPA"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPXA = exp(EV_[1]/pbm.elpar[iel_][1])
        f_   = EXPXA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EXPXA/pbm.elpar[iel_][1]
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = EXPXA/pbm.elpar[iel_][1]/pbm.elpar[iel_][1]
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

