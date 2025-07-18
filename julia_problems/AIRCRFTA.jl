function AIRCRFTA(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : AIRCRFTA
#    *********
# 
#    The aircraft stability problem by Rheinboldt, as a function
#    of the elevator, aileron and rudder deflection controls.
# 
#    Source: Problem 9 in
#    J.J. More',"A collection of nonlinear model problems"
#    Proceedings of the AMS-SIAM Summer Seminar on the Computational
#    Solution of Nonlinear Systems of Equations, Colorado, 1988.
#    Argonne National Laboratory MCS-P60-0289, 1989.
# 
#    SIF input: Ph. Toint, Dec 1989.
# 
#    classification = "C-CNOR2-RN-8-5"
# 
#    Values for the controls
#    1) Elevator
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "AIRCRFTA"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling AIRCRFTA.")
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
        v_["ELVVAL"] = 0.1
        v_["AILVAL"] = 0.0
        v_["RUDVAL"] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        iv,ix_,_ = s2mpj_ii("ROLLRATE",ix_)
        arrset(pb.xnames,iv,"ROLLRATE")
        iv,ix_,_ = s2mpj_ii("PITCHRAT",ix_)
        arrset(pb.xnames,iv,"PITCHRAT")
        iv,ix_,_ = s2mpj_ii("YAWRATE",ix_)
        arrset(pb.xnames,iv,"YAWRATE")
        iv,ix_,_ = s2mpj_ii("ATTCKANG",ix_)
        arrset(pb.xnames,iv,"ATTCKANG")
        iv,ix_,_ = s2mpj_ii("SSLIPANG",ix_)
        arrset(pb.xnames,iv,"SSLIPANG")
        iv,ix_,_ = s2mpj_ii("ELEVATOR",ix_)
        arrset(pb.xnames,iv,"ELEVATOR")
        iv,ix_,_ = s2mpj_ii("AILERON",ix_)
        arrset(pb.xnames,iv,"AILERON")
        iv,ix_,_ = s2mpj_ii("RUDDERDF",ix_)
        arrset(pb.xnames,iv,"RUDDERDF")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("G1",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G1")
        push!(irA,ig)
        push!(icA,ix_["ROLLRATE"])
        push!(valA,Float64(-3.933))
        push!(irA,ig)
        push!(icA,ix_["PITCHRAT"])
        push!(valA,Float64(0.107))
        push!(irA,ig)
        push!(icA,ix_["YAWRATE"])
        push!(valA,Float64(0.126))
        push!(irA,ig)
        push!(icA,ix_["SSLIPANG"])
        push!(valA,Float64(-9.99))
        push!(irA,ig)
        push!(icA,ix_["AILERON"])
        push!(valA,Float64(-45.83))
        push!(irA,ig)
        push!(icA,ix_["RUDDERDF"])
        push!(valA,Float64(-7.64))
        ig,ig_,_ = s2mpj_ii("G2",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G2")
        push!(irA,ig)
        push!(icA,ix_["PITCHRAT"])
        push!(valA,Float64(-0.987))
        push!(irA,ig)
        push!(icA,ix_["ATTCKANG"])
        push!(valA,Float64(-22.95))
        push!(irA,ig)
        push!(icA,ix_["ELEVATOR"])
        push!(valA,Float64(-28.37))
        ig,ig_,_ = s2mpj_ii("G3",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G3")
        push!(irA,ig)
        push!(icA,ix_["ROLLRATE"])
        push!(valA,Float64(0.002))
        push!(irA,ig)
        push!(icA,ix_["YAWRATE"])
        push!(valA,Float64(-0.235))
        push!(irA,ig)
        push!(icA,ix_["SSLIPANG"])
        push!(valA,Float64(5.67))
        push!(irA,ig)
        push!(icA,ix_["AILERON"])
        push!(valA,Float64(-0.921))
        push!(irA,ig)
        push!(icA,ix_["RUDDERDF"])
        push!(valA,Float64(-6.51))
        ig,ig_,_ = s2mpj_ii("G4",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G4")
        push!(irA,ig)
        push!(icA,ix_["PITCHRAT"])
        push!(valA,Float64(1.0))
        push!(irA,ig)
        push!(icA,ix_["ATTCKANG"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["ELEVATOR"])
        push!(valA,Float64(-1.168))
        ig,ig_,_ = s2mpj_ii("G5",ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G5")
        push!(irA,ig)
        push!(icA,ix_["YAWRATE"])
        push!(valA,Float64(-1.0))
        push!(irA,ig)
        push!(icA,ix_["SSLIPANG"])
        push!(valA,Float64(-0.196))
        push!(irA,ig)
        push!(icA,ix_["AILERON"])
        push!(valA,Float64(-0.0071))
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
        pb.xlower[ix_["ELEVATOR"]] = v_["ELVVAL"]
        pb.xupper[ix_["ELEVATOR"]] = v_["ELVVAL"]
        pb.xlower[ix_["AILERON"]] = v_["AILVAL"]
        pb.xupper[ix_["AILERON"]] = v_["AILVAL"]
        pb.xlower[ix_["RUDDERDF"]] = v_["RUDVAL"]
        pb.xupper[ix_["RUDDERDF"]] = v_["RUDVAL"]
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        if haskey(ix_,"ELEVATOR")
            pb.x0[ix_["ELEVATOR"]] = Float64(v_["ELVVAL"])
        else
            pb.y0[findfirst(x->x==ig_["ELEVATOR"],pbm.congrps)] = Float64(v_["ELVVAL"])
        end
        if haskey(ix_,"AILERON")
            pb.x0[ix_["AILERON"]] = Float64(v_["AILVAL"])
        else
            pb.y0[findfirst(x->x==ig_["AILERON"],pbm.congrps)] = Float64(v_["AILVAL"])
        end
        if haskey(ix_,"RUDDERDF")
            pb.x0[ix_["RUDDERDF"]] = Float64(v_["RUDVAL"])
        else
            pb.y0[findfirst(x->x==ig_["RUDDERDF"],pbm.congrps)] = Float64(v_["RUDVAL"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "en2PR", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "E1A"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "PITCHRAT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "YAWRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E1B"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "YAWRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "ATTCKANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E1C"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "ATTCKANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SSLIPANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E1D"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "PITCHRAT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "ATTCKANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2A"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "ROLLRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "YAWRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E2B"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "ROLLRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "SSLIPANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3A"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "ROLLRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "PITCHRAT"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "E3B"
        ie,ie_,newelt = s2mpj_ii(ename,ie_)
        if newelt > 0
            arrset(pbm.elftype,ie,"en2PR")
            arrset(ielftype,ie,iet_["en2PR"])
        end
        vname = "ROLLRATE"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        vname = "ATTCKANG"
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["G1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1A"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.727))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E1B"])
        loaset(pbm.grelw,ig,posel,Float64(8.39))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1C"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-684.4))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E1D"])
        loaset(pbm.grelw,ig,posel,Float64(63.5))
        ig = ig_["G2"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2A"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.949))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E2B"])
        loaset(pbm.grelw,ig,posel,Float64(0.173))
        ig = ig_["G3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3A"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-0.716))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["E3B"])
        loaset(pbm.grelw,ig,posel,Float64(-1.578))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E1D"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.132))
        ig = ig_["G4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E2B"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["G5"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["E3B"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,1.)
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
        pb.pbclass = "C-CNOR2-RN-8-5"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm


    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "en2PR"

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

