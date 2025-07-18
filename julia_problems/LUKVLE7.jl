function LUKVLE7(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : LUKVLE7
#    *********
# 
#    Source: Problem 5.7, A trigonometric tridiagonal function with 
#    simplified five-diagonal constraints, due to L. Luksan and J. Vlcek,
#    "Sparse and partially separable test problems for 
#    unconstrained and equality constrained optimization",
#    Technical Report 767, Inst. Computer Science, Academy of Sciences
#    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
# 
#    SIF input: Nick Gould, April 2001
# 
#    classification = "C-COOR2-AY-V-V"
# 
#    some useful parameters, including N, the number of variables.
# 
#       Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   10000          $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "LUKVLE7"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling LUKVLE7.")
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
        if nargin<1
            v_["N"] = Int64(10);  #  SIF file default value
        else
            v_["N"] = Int64(args[1]);
        end
# IE N                   100000         $-PARAMETER
        v_["1"] = 1
        v_["2"] = 2
        v_["3"] = 3
        v_["4"] = 4
        v_["N-1"] = -1+v_["N"]
        v_["N-2"] = -2+v_["N"]
        v_["N-3"] = -3+v_["N"]
        v_["N"] = Float64(v_["N"])
        v_["-N"] = -1.0*v_["N"]
        v_["N+1"] = 1.0+v_["N"]
        v_["N.N+1"] = v_["N"]*v_["N+1"]
        v_["-N.N+1/2"] = -0.5*v_["N.N+1"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["1"])))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["1"]))])
        push!(valA,Float64(4.0))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["2"]))])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["2"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["2"])))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["2"]))])
        push!(valA,Float64(6.0))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["3"]))])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["3"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["3"])))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["N-1"]))])
        push!(valA,Float64(6.0))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["N-3"]))])
        push!(valA,Float64(-1.0))
        ig,ig_,_ = s2mpj_ii("C"*string(Int64(v_["4"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"C"*string(Int64(v_["4"])))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["N"]))])
        push!(valA,Float64(2.0))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["N-2"]))])
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
        pbm.gconst[ig_["OBJ"]] = Float64(v_["-N.N+1/2"])
        pbm.gconst[ig_["C"*string(Int64(v_["2"]))]] = Float64(2.0)
        pbm.gconst[ig_["C"*string(Int64(v_["3"]))]] = Float64(2.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["X"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSIN", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eCOS", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eCUBEP", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "SI"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSIN")
            arrset(ielftype,ie,iet_["eSIN"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "CO"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOS")
            arrset(ielftype,ie,iet_["eCOS"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "C1"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C1"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C2"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C2"*string(Int64(v_["1"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCUBEP")
        arrset(ielftype,ie,iet_["eCUBEP"])
        ename = "C1"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C2"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C2"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["3"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C3"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C3"*string(Int64(v_["2"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["4"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCUBEP")
        arrset(ielftype,ie,iet_["eCUBEP"])
        ename = "C1"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C2"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C2"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C3"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C3"*string(Int64(v_["3"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N-2"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eCUBEP")
        arrset(ielftype,ie,iet_["eCUBEP"])
        ename = "C1"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C1"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="W",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "C2"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        ename = "C2"*string(Int64(v_["4"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "X"*string(Int64(v_["N-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="V",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CO"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["SI"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            v_["I-1"] = -1+I
            v_["I"] = Float64(I)
            v_["-I"] = -1.0*I
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CO"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-I"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["SI"*string(Int64(v_["I-1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(I))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["SI"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-I"]))
        end
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["CO"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["-N"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["SI"*string(Int64(v_["N-1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["N"]))
        ig = ig_["C"*string(Int64(v_["1"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-4.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2"*string(Int64(v_["1"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C"*string(Int64(v_["2"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(8.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-4.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C3"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["C"*string(Int64(v_["3"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(8.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-4.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C3"*string(Int64(v_["3"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        ig = ig_["C"*string(Int64(v_["4"]))]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C1"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(8.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C2"*string(Int64(v_["4"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               -2.2421E+02
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COOR2-AY-V-V"
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

    elseif action == "eSIN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SINV = sin(EV_[1])
        f_   = SINV
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = cos(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -SINV
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eCOS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        COSV = cos(EV_[1])
        f_   = COSV
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -sin(EV_[1])
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -COSV
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]*EV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]
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

    elseif action == "eCUBEP"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = EV_[1]^3-EV_[1]*EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 3.0*EV_[1]^2-EV_[2]
            g_[2] = -EV_[1]
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 6.0*EV_[1]
                H_[1,2] = -1.0
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

