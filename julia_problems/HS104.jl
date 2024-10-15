function HS104(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS104
#    *********
# 
#    Source: problem 104 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, August 1991.
# 
#    classification = "C-OOR2-AN-8-5"
# 
#    Number of variables
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS104"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 8
        v_["1"] = 1
        v_["4"] = 4
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0e+0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0e+0)
        ig,ig_,_ = s2mpj_ii("C1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C1")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0e-1)
        ig,ig_,_ = s2mpj_ii("C2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C2")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(1.0e-1)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(1.0e-1)
        ig,ig_,_ = s2mpj_ii("C3",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C3")
        ig,ig_,_ = s2mpj_ii("C4",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"C4")
        ig,ig_,_ = s2mpj_ii("C5",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"C5")
        iv = ix_["X1"]
        pbm.A[ig,iv] += Float64(-1.0e+0)
        iv = ix_["X2"]
        pbm.A[ig,iv] += Float64(-1.0e+0)
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
        #%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = fill(1.0e+0,ngrp)
        pbm.gconst[ig_["OBJ"]] = Float64(-1.0e+1)
        pbm.gconst[ig_["C5"]] = Float64(-9.0e+0)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        grange[gegrps,1] = fill(Inf,pb.nge)
        arrset(grange,ig_["C5"],Float64(3.2e+0))
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = fill(1.0e-1,pb.n)
        pb.xupper = fill(1.0e+1,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(3.0)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(3.0)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(0.4)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(0.4)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(0.2)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(0.2)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"X6")
            pb.x0[ix_["X6"]] = Float64(6.0)
        else
            pb.y0[findfirst(x->x==ig_["X6"],pbm.congrps)] = Float64(6.0)
        end
        if haskey(ix_,"X7")
            pb.x0[ix_["X7"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["X7"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"X8")
            pb.x0[ix_["X8"]] = Float64(0.5)
        else
            pb.y0[findfirst(x->x==ig_["X8"],pbm.congrps)] = Float64(0.5)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"P1")
        loaset(elftp,it,2,"P2")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "OE1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X1"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.67))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.67))
        ename = "OE2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X2"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(0.67))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.67))
        ename = "C1E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "C2E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "C3E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        ename = "C3E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X5"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.71))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        ename = "C3E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X3"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X7"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.3))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        ename = "C4E1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        ename = "C4E2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X6"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-0.71))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.0))
        ename = "C4E3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        vname = "X4"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "X8"
        iv,ix_,pb  = (
              s2mpj_nlx(vname,ix_,pb,1,Float64(1.0e-1),Float64(1.0e+1),nothing))
        posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        posep = findfirst(x->x=="P1",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(-1.3))
        posep = findfirst(x->x=="P2",elftp[ielftype[ie]])
        loaset(pbm.elpar,ie,posep,Float64(1.0))
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.0e-1))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE2"])
        loaset(pbm.grelw,ig,posel,Float64(4.0e-1))
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.88e-2))
        ig = ig_["C2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.88e-2))
        ig = ig_["C3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C3E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.0e+0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C3E2"])
        loaset(pbm.grelw,ig,posel,Float64(2.0e+0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C3E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.88e-2))
        ig = ig_["C4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C4E1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.0e+0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C4E2"])
        loaset(pbm.grelw,ig,posel,Float64(2.0e+0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C4E3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(5.88e-2))
        ig = ig_["C5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["OE1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(4.0e-1))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["OE2"])
        loaset(pbm.grelw,ig,posel,Float64(4.0e-1))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               3.9511634396
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[1:pb.nle] = grange[legrps]
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-OOR2-AN-8-5"
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

    elseif action == "ePROD"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^pbm.elpar[iel_][2])
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1]  = (
                  pbm.elpar[iel_][1]*(EV_[1]^(pbm.elpar[iel_][1]-1.0))*(EV_[2]^pbm.elpar[iel_][2]))
            g_[2]  = (
                  pbm.elpar[iel_][2]*(EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^(pbm.elpar[iel_][2]-1.0)))
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1]  = (
                      pbm.elpar[iel_][1]*(EV_[1]^(pbm.elpar[iel_][1]-2.0))*(pbm.elpar[iel_][1]-1.0)*(EV_[2]^pbm.elpar[iel_][2]))
                H_[1,2]  = (
                      pbm.elpar[iel_][1]*(EV_[1]^(pbm.elpar[iel_][1]-1.0))*pbm.elpar[iel_][2]*(EV_[2]^(pbm.elpar[iel_][2]-1.0)))
                H_[2,1] = H_[1,2]
                H_[2,2]  = (
                      pbm.elpar[iel_][2]*(pbm.elpar[iel_][2]-1.0)*(EV_[1]^pbm.elpar[iel_][1])*(EV_[2]^(pbm.elpar[iel_][2]-2.0)))
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

