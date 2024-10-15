function POLAK2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : POLAK2
#    *********
# 
#    A nonlinear minmax problem in ten variables.
# 
#    Source: 
#    E. Polak, D.H. Mayne and J.E. Higgins,
#    "Superlinearly convergent algorithm for min-max problems"
#    JOTA 69, pp. 407-439, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-LOR2-AN-11-2"
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "POLAK2"

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
        v_["10"] = 10
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["10"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
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
        ig,ig_,_ = s2mpj_ii("F2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"F2")
        iv = ix_["U"]
        pbm.A[ig,iv] += Float64(-1.0)
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
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.1),pb.n)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(100.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(100.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eEL", iet_)
        loaset(elftv,it,1,"XX1")
        loaset(elftv,it,2,"XX2")
        loaset(elftv,it,3,"XX3")
        loaset(elftv,it,4,"XX4")
        loaset(elftv,it,5,"XX5")
        loaset(elftv,it,6,"XX6")
        loaset(elftv,it,7,"XX7")
        loaset(elftv,it,8,"XX8")
        loaset(elftv,it,9,"XX9")
        loaset(elftv,it,10,"XX10")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL")
        arrset(ielftype,ie,iet_["eEL"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX8",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX9",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX10",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(2.0))
        ename = "E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eEL")
        arrset(ielftype,ie,iet_["eEL"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX3",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX4",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX5",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX6",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX7",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX8",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX9",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.1))
        posev = findfirst(x->x=="XX10",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-2.0))
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
        ig = ig_["F2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               54.598146
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
        pb.pbclass = "C-LOR2-AN-11-2"
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

    elseif action == "eEL"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        A = 1.0e-8*EV_[1]*EV_[1]+(EV_[2]+pbm.elpar[iel_][1])^2
        A = A+EV_[3]*EV_[3]+4.0*EV_[4]*EV_[4]
        A = A+EV_[5]*EV_[5]+EV_[6]*EV_[6]+EV_[7]*EV_[7]
        A = A+EV_[8]*EV_[8]+EV_[9]*EV_[9]+EV_[10]*EV_[10]
        EA = exp(A)
        f_   = EA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0e-8*EV_[1]*EA
            g_[2] = 2.0*(EV_[2]+pbm.elpar[iel_][1])*EA
            g_[3] = 2.0*EV_[3]*EA
            g_[4] = 8.0*EV_[4]*EA
            g_[5] = 2.0*EV_[5]*EA
            g_[6] = 2.0*EV_[6]*EA
            g_[7] = 2.0*EV_[7]*EA
            g_[8] = 2.0*EV_[8]*EA
            g_[9] = 2.0*EV_[9]*EA
            g_[10] = 2.0*EV_[10]*EA
            if nargout>2
                H_ = zeros(Float64,10,10)
                H_[1,1] = 2.0e-8*EA*(1.0+2.0e-8*EV_[1]^2)
                H_[1,2] = 4.0e-8*EV_[1]*(EV_[2]+pbm.elpar[iel_][1])*EA
                H_[2,1] = H_[1,2]
                H_[1,3] = 4.0e-8*EV_[1]*EV_[3]*EA
                H_[3,1] = H_[1,3]
                H_[1,4] = 1.6e-7*EV_[1]*EV_[4]*EA
                H_[4,1] = H_[1,4]
                H_[1,5] = 4.0e-8*EV_[1]*EV_[5]*EA
                H_[5,1] = H_[1,5]
                H_[1,6] = 4.0e-8*EV_[1]*EV_[6]*EA
                H_[6,1] = H_[1,6]
                H_[1,7] = 4.0e-8*EV_[1]*EV_[7]*EA
                H_[7,1] = H_[1,7]
                H_[1,8] = 4.0e-8*EV_[1]*EV_[8]*EA
                H_[8,1] = H_[1,8]
                H_[1,9] = 4.0e-8*EV_[1]*EV_[9]*EA
                H_[9,1] = H_[1,9]
                H_[1,10] = 4.0e-8*EV_[1]*EV_[10]*EA
                H_[10,1] = H_[1,10]
                H_[2,2] = 2.0*EA*(1.0+2.0*(EV_[2]+pbm.elpar[iel_][1])^2)
                H_[2,3] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[3]*EA
                H_[3,2] = H_[2,3]
                H_[2,4] = 16.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[4]*EA
                H_[4,2] = H_[2,4]
                H_[2,5] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[5]*EA
                H_[5,2] = H_[2,5]
                H_[2,6] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[6]*EA
                H_[6,2] = H_[2,6]
                H_[2,7] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[7]*EA
                H_[7,2] = H_[2,7]
                H_[2,8] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[8]*EA
                H_[8,2] = H_[2,8]
                H_[2,9] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[9]*EA
                H_[9,2] = H_[2,9]
                H_[2,10] = 4.0*(EV_[2]+pbm.elpar[iel_][1])*EV_[10]*EA
                H_[10,2] = H_[2,10]
                H_[3,3] = 2.0*EA*(1.0+2.0*EV_[3]*EV_[3])
                H_[3,4] = 16.0*EV_[3]*EV_[4]*EA
                H_[4,3] = H_[3,4]
                H_[3,5] = 4.0*EV_[3]*EV_[5]*EA
                H_[5,3] = H_[3,5]
                H_[3,6] = 4.0*EV_[3]*EV_[6]*EA
                H_[6,3] = H_[3,6]
                H_[3,7] = 4.0*EV_[3]*EV_[7]*EA
                H_[7,3] = H_[3,7]
                H_[3,8] = 4.0*EV_[3]*EV_[8]*EA
                H_[8,3] = H_[3,8]
                H_[3,9] = 4.0*EV_[3]*EV_[9]*EA
                H_[9,3] = H_[3,9]
                H_[3,10] = 4.0*EV_[3]*EV_[10]*EA
                H_[10,3] = H_[3,10]
                H_[4,4] = 8.0*EA*(1.0+8.0*EV_[4]*EV_[4])
                H_[4,5] = 16.0*EV_[4]*EV_[5]*EA
                H_[5,4] = H_[4,5]
                H_[4,6] = 16.0*EV_[4]*EV_[6]*EA
                H_[6,4] = H_[4,6]
                H_[4,7] = 16.0*EV_[4]*EV_[7]*EA
                H_[7,4] = H_[4,7]
                H_[4,8] = 16.0*EV_[4]*EV_[8]*EA
                H_[8,4] = H_[4,8]
                H_[4,9] = 16.0*EV_[4]*EV_[9]*EA
                H_[9,4] = H_[4,9]
                H_[4,10] = 16.0*EV_[4]*EV_[10]*EA
                H_[10,4] = H_[4,10]
                H_[5,5] = 2.0*EA*(1.0+2.0*EV_[5]*EV_[5])
                H_[5,6] = 4.0*EV_[5]*EV_[6]*EA
                H_[6,5] = H_[5,6]
                H_[5,7] = 4.0*EV_[5]*EV_[7]*EA
                H_[7,5] = H_[5,7]
                H_[5,8] = 4.0*EV_[5]*EV_[8]*EA
                H_[8,5] = H_[5,8]
                H_[5,9] = 4.0*EV_[5]*EV_[9]*EA
                H_[9,5] = H_[5,9]
                H_[5,10] = 4.0*EV_[5]*EV_[10]*EA
                H_[10,5] = H_[5,10]
                H_[6,6] = 2.0*EA*(1.0+2.0*EV_[6]*EV_[6])
                H_[6,7] = 4.0*EV_[6]*EV_[7]*EA
                H_[7,6] = H_[6,7]
                H_[6,8] = 4.0*EV_[6]*EV_[8]*EA
                H_[8,6] = H_[6,8]
                H_[6,9] = 4.0*EV_[6]*EV_[9]*EA
                H_[9,6] = H_[6,9]
                H_[6,10] = 4.0*EV_[6]*EV_[10]*EA
                H_[10,6] = H_[6,10]
                H_[7,7] = 2.0*EA*(1.0+2.0*EV_[7]*EV_[7])
                H_[7,8] = 4.0*EV_[7]*EV_[8]*EA
                H_[8,7] = H_[7,8]
                H_[7,9] = 4.0*EV_[7]*EV_[9]*EA
                H_[9,7] = H_[7,9]
                H_[7,10] = 4.0*EV_[7]*EV_[10]*EA
                H_[10,7] = H_[7,10]
                H_[8,8] = 2.0*EA*(1.0+2.0*EV_[8]*EV_[8])
                H_[8,9] = 4.0*EV_[8]*EV_[9]*EA
                H_[9,8] = H_[8,9]
                H_[8,10] = 4.0*EV_[8]*EV_[10]*EA
                H_[10,8] = H_[8,10]
                H_[9,9] = 2.0*EA*(1.0+2.0*EV_[9]*EV_[9])
                H_[9,10] = 4.0*EV_[9]*EV_[10]*EA
                H_[10,9] = H_[9,10]
                H_[10,10] = 2.0*EA*(1.0+2.0*EV_[10]*EV_[10])
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

