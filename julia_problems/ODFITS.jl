function ODFITS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    A simple Origin/Destination matrix fit using a minimum entropy
#    approach.  The objective is a combination of different aims, namely
#    to be close to an a priori matrix for some entries, to be consistent
#    with some traffic counts (for some entries) and to be small (for entries
#    where nothing else is known).
# 
#    The objective function is of the form
#         SUM   m T [ ln( T / a ) - 1 ] + E   SUM  T [ ln ( T  ) - 1 ]
#        i in I  i i       i   i            i in J  i        i
#                +  g   SUM   q  F [ ln( F / c ) - 1 ]
#                     i in K   i  i       i   i
#    with the constraints that all Ti and Fi be positive and that
#                         F  =  SUM p   T
#                          i     j   ij  j
#    where the pij represent path weights from an a priori assignment.
# 
#    Source: a modification of an example in
#    L.G. Willumsen,
#    "Origin-Destination Matrix: static estimation"
#    in "Concise Encyclopedia of Traffic and Transportation Systems"
#    (M. Papageorgiou, ed.), Pergamon Press, 1991.
# 
#    M. Bierlaire, private communication, 1991.
# 
#    SIF input: Ph Toint, Dec 1991.
# 
#    classification = "C-OLR2-MN-10-6"
# 
#    Number of available traffic counts
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ODFITS"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["ARCS"] = 6
        v_["TC1"] = 100.0
        v_["TC2"] = 500.0
        v_["TC3"] = 400.0
        v_["TC4"] = 1100.0
        v_["TC5"] = 600.0
        v_["TC6"] = 700.0
        v_["QLT1"] = 1.0
        v_["QLT2"] = 1.0
        v_["QLT3"] = 1.0
        v_["QLT4"] = 1.0
        v_["QLT5"] = 1.0
        v_["QLT6"] = 1.0
        v_["P131"] = 1.0
        v_["P132"] = 0.0
        v_["P133"] = 0.0
        v_["P134"] = 0.0
        v_["P135"] = 0.0
        v_["P136"] = 0.0
        v_["P141"] = 0.0
        v_["P142"] = 1.0
        v_["P143"] = 0.0
        v_["P144"] = 1.0
        v_["P145"] = 0.0
        v_["P146"] = 0.0
        v_["P231"] = 0.0
        v_["P232"] = 0.0
        v_["P233"] = 1.0
        v_["P234"] = 1.0
        v_["P235"] = 1.0
        v_["P236"] = 0.0
        v_["P241"] = 0.0
        v_["P242"] = 0.0
        v_["P243"] = 0.0
        v_["P244"] = 1.0
        v_["P245"] = 1.0
        v_["P246"] = 1.0
        v_["APV13"] = 90.0
        v_["APV14"] = 450.0
        v_["APV23"] = 360.0
        v_["MU13"] = 0.5
        v_["MU14"] = 0.5
        v_["MU23"] = 0.5
        v_["1/MU13"] = 1.0/v_["MU13"]
        v_["1/MU14"] = 1.0/v_["MU14"]
        v_["1/MU23"] = 1.0/v_["MU23"]
        v_["GAMMA"] = 1.5
        v_["ENTROP"] = 0.2
        v_["1/ENTR"] = 1.0/v_["ENTROP"]
        v_["1"] = 1
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            v_["1/QLT"*string(I)] = 1.0/v_["QLT"*string(I)]
            v_["G/QLT"*string(I)] = v_["1/QLT"*string(I)]*v_["GAMMA"]
        end
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("T13",ix_)
        arrset(pb.xnames,iv,"T13")
        iv,ix_,_ = s2mpj_ii("T14",ix_)
        arrset(pb.xnames,iv,"T14")
        iv,ix_,_ = s2mpj_ii("T23",ix_)
        arrset(pb.xnames,iv,"T23")
        iv,ix_,_ = s2mpj_ii("T24",ix_)
        arrset(pb.xnames,iv,"T24")
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            iv,ix_,_ = s2mpj_ii("F"*string(I),ix_)
            arrset(pb.xnames,iv,"F"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("AP13",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T13"]
        pbm.A[ig,iv] += Float64(-1.0)
        arrset(pbm.gscale,ig,Float64(v_["1/MU13"]))
        ig,ig_,_ = s2mpj_ii("AP14",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T14"]
        pbm.A[ig,iv] += Float64(-1.0)
        arrset(pbm.gscale,ig,Float64(v_["1/MU14"]))
        ig,ig_,_ = s2mpj_ii("AP23",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T23"]
        pbm.A[ig,iv] += Float64(-1.0)
        arrset(pbm.gscale,ig,Float64(v_["1/MU23"]))
        ig,ig_,_ = s2mpj_ii("AP24",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["T24"]
        pbm.A[ig,iv] += Float64(-1.0)
        arrset(pbm.gscale,ig,Float64(v_["1/ENTR"]))
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            ig,ig_,_ = s2mpj_ii("CP"*string(I),ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["F"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            arrset(pbm.gscale,ig,Float64(v_["G/QLT"*string(I)]))
        end
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(I))
            iv = ix_["F"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["T13"]
            pbm.A[ig,iv] += Float64(v_["P13"*string(I)])
            iv = ix_["T14"]
            pbm.A[ig,iv] += Float64(v_["P14"*string(I)])
            iv = ix_["T23"]
            pbm.A[ig,iv] += Float64(v_["P23"*string(I)])
            iv = ix_["T24"]
            pbm.A[ig,iv] += Float64(v_["P24"*string(I)])
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(0.1,pb.n)
        pb.xupper = fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["T13"]] = Float64(v_["APV13"])
        pb.x0[ix_["T14"]] = Float64(v_["APV14"])
        pb.x0[ix_["T23"]] = Float64(v_["APV23"])
        pb.x0[ix_["T24"]] = Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            pb.x0[ix_["F"*string(I)]] = Float64(v_["TC"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXLOGX", iet_)
        loaset(elftv,it,1,"X")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"DEN")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "TFIT13"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXLOGX")
        arrset(ielftype,ie,iet_["eXLOGX"])
        vname = "T13"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="DEN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["APV13"]))
        ename = "TFIT23"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXLOGX")
        arrset(ielftype,ie,iet_["eXLOGX"])
        vname = "T23"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="DEN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["APV23"]))
        ename = "TFIT14"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXLOGX")
        arrset(ielftype,ie,iet_["eXLOGX"])
        vname = "T14"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="DEN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(v_["APV14"]))
        ename = "TFIT24"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXLOGX")
        arrset(ielftype,ie,iet_["eXLOGX"])
        vname = "T24"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="DEN",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            ename = "CFIT"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eXLOGX")
            arrset(ielftype,ie,iet_["eXLOGX"])
            vname = "F"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,Float64(0.1),nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="DEN",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["TC"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["AP13"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TFIT13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AP14"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TFIT14"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AP23"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TFIT23"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["AP24"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["TFIT24"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        for I = Int64(v_["1"]):Int64(v_["ARCS"])
            ig = ig_["CP"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CFIT"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO ODFITS             -2380.026775
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
        pb.pbclass = "C-OLR2-MN-10-6"
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

    elseif action == "eXLOGX"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        LOGX = log(EV_[1]/pbm.elpar[iel_][1])
        f_   = EV_[1]*LOGX
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0+LOGX
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 1.0/EV_[1]
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

