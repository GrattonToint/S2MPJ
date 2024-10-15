function HONG(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    Source: Se June Hong/Chid Apte
# 
#    SIF input: A.R.Conn, Jan 1991.
# 
#    classification = "C-OLR2-AN-4-1"
# 
#   Problem parameters
# 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HONG"

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
        iv,ix_,_ = s2mpj_ii("T1",ix_)
        arrset(pb.xnames,iv,"T1")
        iv,ix_,_ = s2mpj_ii("T2",ix_)
        arrset(pb.xnames,iv,"T2")
        iv,ix_,_ = s2mpj_ii("T3",ix_)
        arrset(pb.xnames,iv,"T3")
        iv,ix_,_ = s2mpj_ii("T4",ix_)
        arrset(pb.xnames,iv,"T4")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("SUM1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"SUM1")
        iv = ix_["T1"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["T2"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["T3"]
        pbm.A[ig,iv] += Float64(1.0)
        iv = ix_["T4"]
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
        pbm.gconst[ig_["SUM1"]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.0,pb.n)
        pb.xupper = fill(1.0,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEXP", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        loaset(elftp,it,3,"P3")
        loaset(elftp,it,4,"P4")
        loaset(elftp,it,5,"P5")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
        end
        vname = "T1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(25.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.92))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.08))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.38))
        ename = "E2"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
        end
        vname = "T2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(50.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.95))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(3.95))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.11))
        ename = "E3"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
        end
        vname = "T3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(9.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-4.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.66))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1657834.0))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.48))
        ename = "E4"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eEXP")
            arrset(ielftype,ie,iet_["eEXP"])
        end
        vname = "T4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.0),Float64(1.0),Float64(0.5))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(20000.0))
        posep = findfirst(x->x=="P3",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.11))
        posep = findfirst(x->x=="P4",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.89))
        posep = findfirst(x->x=="P5",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.00035))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E4"])
        loaset(pbm.grelw,ig,posel, 1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution unknown
        pb.objlower = -4.0
        pb.objupper = 300.0
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
        pb.pbclass = "C-OLR2-AN-4-1"
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

    elseif action == "eEXP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        XTOT = pbm.elpar[iel_][1]+pbm.elpar[iel_][2]*EV_[1]
        EP5 = exp(pbm.elpar[iel_][5]*XTOT)
        f_   = pbm.elpar[iel_][3]+pbm.elpar[iel_][4]*EP5
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = pbm.elpar[iel_][2]*pbm.elpar[iel_][4]*pbm.elpar[iel_][5]*EP5
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1]  = (
                      pbm.elpar[iel_][2]*pbm.elpar[iel_][2]*pbm.elpar[iel_][4]*pbm.elpar[iel_][5]*pbm.elpar[iel_][5]*EP5)
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

