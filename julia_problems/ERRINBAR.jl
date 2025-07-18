function ERRINBAR(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ERRINBAR
#    *********
#    A nonlinearly constrained problem.
# 
#    Source:
#    Ph. Toint,
#    an error in specifying the proroblem TENBARS1
#    Oct. 90.
# 
#    SIF input: Ph. Toint, October 1990.
# 
#    classification = "C-CLOR2-MY-18-9"
# 
#    Constants 
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ERRINBAR"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling ERRINBAR.")
    end

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
        v_["2"] = 2
        v_["8"] = 8
        v_["10"] = 10
        v_["SQ2"] = sqrt(2.0)
        v_["SQ8"] = sqrt(8.0)
        v_["1/SQ8"] = 1.0/v_["SQ8"]
        v_["C0"] = 2.53106
        v_["-P0"] = -589.884
        v_["C0SQ2"] = v_["C0"]*v_["SQ2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["8"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["10"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        push!(irA,ig)
        push!(icA,ix_["X1"])
        push!(valA,Float64(v_["C0"]))
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(v_["C0SQ2"]))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(v_["C0SQ2"]))
        push!(irA,ig)
        push!(icA,ix_["X4"])
        push!(valA,Float64(v_["C0"]))
        push!(irA,ig)
        push!(icA,ix_["X5"])
        push!(valA,Float64(v_["C0"]))
        push!(irA,ig)
        push!(icA,ix_["X6"])
        push!(valA,Float64(v_["C0"]))
        push!(irA,ig)
        push!(icA,ix_["X7"])
        push!(valA,Float64(v_["C0SQ2"]))
        push!(irA,ig)
        push!(icA,ix_["X8"])
        push!(valA,Float64(v_["C0SQ2"]))
        push!(irA,ig)
        push!(icA,ix_["X9"])
        push!(valA,Float64(v_["C0"]))
        push!(irA,ig)
        push!(icA,ix_["X10"])
        push!(valA,Float64(v_["C0"]))
        for I = Int64(v_["1"]):Int64(v_["8"])
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"C"*string(I))
        end
        ig,ig_,_ = s2mpj_ii("STR",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"STR")
        push!(irA,ig)
        push!(icA,ix_["U2"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["U4"])
        push!(valA,Float64(-1.0))
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
        pbm.gconst[ig_["C4"]] = Float64(v_["-P0"])
        pbm.gconst[ig_["C8"]] = Float64(v_["-P0"])
        pbm.gconst[ig_["STR"]] = Float64(22.86)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["10"])
            pb.xlower[ix_["X"*string(I)]] = 0.645
        end
        for I = Int64(v_["2"]):Int64(v_["2"]):Int64(v_["8"])
            pb.xlower[ix_["U"*string(I)]] = -50.8
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eXU", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"U")
        it,iet_,_ = s2mpj_ii( "eXUPV", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"U")
        loaset(elftv,it,3,"V")
        it,iet_,_ = s2mpj_ii( "eXUMV", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"U")
        loaset(elftv,it,3,"V")
        it,iet_,_ = s2mpj_ii( "eXBIG", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"U")
        loaset(elftv,it,3,"V")
        loaset(elftv,it,4,"W")
        loaset(elftv,it,5,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "EA"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXU")
        arrset(ielftype,ie,iet_["eXU"])
        vname = "X1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EB"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXUPV")
        arrset(ielftype,ie,iet_["eXUPV"])
        vname = "X3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EC"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eXUMV")
            arrset(ielftype,ie,iet_["eXUMV"])
        end
        vname = "X6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "ED"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXBIG")
        arrset(ielftype,ie,iet_["eXBIG"])
        vname = "X7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U1"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EE"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eXUMV")
            arrset(ielftype,ie,iet_["eXUMV"])
        end
        vname = "X5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EF"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eXUMV")
            arrset(ielftype,ie,iet_["eXUMV"])
        end
        vname = "X2"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EG"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXU")
        arrset(ielftype,ie,iet_["eXU"])
        vname = "X4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EH"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eXBIG")
        arrset(ielftype,ie,iet_["eXBIG"])
        vname = "X8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U4"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U5"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EI"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eXUMV")
            arrset(ielftype,ie,iet_["eXUMV"])
        end
        vname = "X9"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U3"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U7"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "EJ"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"eXUMV")
            arrset(ielftype,ie,iet_["eXUMV"])
        end
        vname = "X10"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U6"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="U",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "U8"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["C1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EA"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EC"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EB"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ED"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        ig = ig_["C2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EB"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["EE"])
        loaset(pbm.grelw,ig,posel, 1.)
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ED"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        ig = ig_["C3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EF"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EG"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EI"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EF"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EE"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EC"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C6"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EH"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EJ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C7"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ED"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EI"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        ig = ig_["C8"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["ED"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["1/SQ8"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["EJ"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               28.04525
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-CLOR2-MY-18-9"
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

    elseif action == "eXU"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = EV_[2]
            g_[2] = EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = 1.0
                H_[2,1] = H_[1,2]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eXUPV"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[2,1] = U_[2,1]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[2]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = IV_[1]
            g_[1] = IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = 1.0
                H_[1,2] = H_[2,1]
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

    elseif action == "eXUMV"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,3)
        IV_ =  zeros(Float64,2)
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]-1
        U_[2,1] = U_[2,1]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[2]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = IV_[1]
            g_[1] = IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = 1.0
                H_[1,2] = H_[2,1]
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

    elseif action == "eXBIG"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,2,5)
        IV_ =  zeros(Float64,2)
        U_[1,2] = U_[1,2]+1
        U_[1,3] = U_[1,3]+1
        U_[1,4] = U_[1,4]-1
        U_[1,5] = U_[1,5]-1
        U_[2,1] = U_[2,1]+1
        IV_[1] = dot(U_[1,:],EV_)
        IV_[2] = dot(U_[2,:],EV_)
        f_   = IV_[2]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[2] = IV_[1]
            g_[1] = IV_[2]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[2,1] = 1.0
                H_[1,2] = H_[2,1]
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

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

