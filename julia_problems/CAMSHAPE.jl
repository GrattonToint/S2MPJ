function CAMSHAPE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CAMSHAPE
#    *********
# 
#    Maximize the area of the valve opening for one rotation of a convex cam 
#    with constraints on the curvature and on the radius of the cam
# 
#    This is problem 4 in the COPS (Version 2) collection of 
#    E. Dolan and J. More'
#    see "Benchmarking Optimization Software with COPS"
#    Argonne National Labs Technical Report ANL/MCS-246 (2000)
# 
#    SIF input: Nick Gould, November 2000
# 
#    classification = "C-LOR2-AN-V-V"
# 
#    The number of discretization points
# 
#       Alternative values for the SIF file parameters:
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# IE N                   400            $-PARAMETER
# IE N                   800            $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CAMSHAPE"

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
        v_["RV"] = 1.0
        v_["RMAX"] = 2.0
        v_["RMIN"] = 1.0
        v_["RAV"] = v_["RMIN"]+v_["RMAX"]
        v_["RAV"] = 0.5*v_["RAV"]
        v_["PI/4"] = atan(1.0)
        v_["PI"] = 4.0*v_["PI/4"]
        v_["ALPHA"] = 1.5
        v_["N+1"] = 1+v_["N"]
        v_["5(N+1)"] = 5*v_["N+1"]
        v_["5(N+1)"] = Float64(v_["5(N+1)"])
        v_["DTHETA"] = 2.0*v_["PI"]
        v_["DTHETA"] = v_["DTHETA"]/v_["5(N+1)"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["PIRV"] = v_["PI"]*v_["RV"]
        v_["PIRV/N"] = v_["PIRV"]/v_["RN"]
        v_["-PIRV/N"] = -1.0*v_["PIRV/N"]
        v_["CDTHETA"] = cos(v_["DTHETA"])
        v_["2CDTHETA"] = 2.0*v_["CDTHETA"]
        v_["ADTHETA"] = v_["ALPHA"]*v_["DTHETA"]
        v_["-ADTHETA"] = -1.0*v_["ADTHETA"]
        v_["2ADTHETA"] = 2.0*v_["ADTHETA"]
        v_["-RMIN"] = -1.0*v_["RMIN"]
        v_["-RMAX"] = -1.0*v_["RMAX"]
        v_["-2RMAX"] = -2.0*v_["RMAX"]
        v_["RMIN2"] = v_["RMIN"]*v_["RMIN"]
        v_["RMIN2CD"] = v_["RMIN"]*v_["2CDTHETA"]
        v_["RMAX2CD"] = v_["RMAX"]*v_["2CDTHETA"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("R"*string(I),ix_)
            arrset(pb.xnames,iv,"R"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("AREA",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-PIRV/N"])
        end
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            ig,ig_,_ = s2mpj_ii("CO"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"CO"*string(I))
        end
        ig,ig_,_ = s2mpj_ii("E1",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E1")
        iv = ix_["R"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["-RMIN"])
        iv = ix_["R"*string(Int64(v_["2"]))]
        pbm.A[ig,iv] += Float64(v_["RMIN2CD"])
        v_["R"] = v_["RMIN2CD"]-v_["RMIN"]
        ig,ig_,_ = s2mpj_ii("E2",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E2")
        iv = ix_["R"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(v_["R"])
        ig,ig_,_ = s2mpj_ii("E3",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E3")
        iv = ix_["R"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["-RMAX"])
        iv = ix_["R"*string(Int64(v_["N-1"]))]
        pbm.A[ig,iv] += Float64(v_["RMAX2CD"])
        ig,ig_,_ = s2mpj_ii("E4",ig_)
        arrset(gtype,ig,"<=")
        arrset(pb.cnames,ig,"E4")
        iv = ix_["R"*string(Int64(v_["N"]))]
        pbm.A[ig,iv] += Float64(v_["-2RMAX"])
        ig,ig_,_ = s2mpj_ii("CU"*string(Int64(v_["0"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CU"*string(Int64(v_["0"])))
        iv = ix_["R"*string(Int64(v_["1"]))]
        pbm.A[ig,iv] += Float64(1.0)
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("CU"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CU"*string(I))
            iv = ix_["R"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["R"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
        end
        ig,ig_,_ = s2mpj_ii("CU"*string(Int64(v_["N"])),ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CU"*string(Int64(v_["N"])))
        iv = ix_["R"*string(Int64(v_["N"]))]
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        pbm.gconst[ig_["E2"]] = Float64(v_["RMIN2"])
        v_["R"] = v_["-ADTHETA"]+v_["RMIN"]
        pbm.gconst[ig_["CU"*string(Int64(v_["0"]))]] = Float64(v_["R"])
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pbm.gconst[ig_["CU"*string(I)]] = Float64(v_["-ADTHETA"])
        end
        v_["R"] = v_["-ADTHETA"]-v_["RMAX"]
        pbm.gconst[ig_["CU"*string(Int64(v_["N"]))]] = Float64(v_["R"])
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[legrps,1] = fill(Inf,pb.nle)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for I = Int64(v_["0"]):Int64(v_["N"])
            arrset(grange,ig_["CU"*string(I)],Float64(v_["2ADTHETA"]))
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["R"*string(I)]] = v_["RMIN"]
            pb.xupper[ix_["R"*string(I)]] = v_["RMAX"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.x0[ix_["R"*string(I)]] = Float64(v_["RAV"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQR", iet_)
        loaset(elftv,it,1,"X")
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I-1"] = -1+I
            v_["I+1"] = 1+I
            ename = "RA"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "R"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "R"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "RB"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "R"*string(Int64(v_["I+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "R"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        ename = "RA"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"ePROD")
        arrset(ielftype,ie,iet_["ePROD"])
        ename = "RA"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "R"*string(Int64(v_["N"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "RA"*string(Int64(v_["N"]))
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        vname = "R"*string(Int64(v_["N-1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "R2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQR")
        arrset(ielftype,ie,iet_["eSQR"])
        vname = "R"*string(Int64(v_["N"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig = ig_["CO"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["RA"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["RA"*string(Int64(v_["I+1"]))])
            loaset(pbm.grelw,ig,posel,Float64(-1.0))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["RB"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["2CDTHETA"]))
        end
        ig = ig_["E1"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["RA"*string(Int64(v_["2"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["E3"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["RA"*string(Int64(v_["N"]))])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        ig = ig_["E4"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["R2"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["2CDTHETA"]))
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLUTION             -4.2841D+00   $ (NH=100)
# LO SOLUTION             -4.2785D+00   $ (NH=200)
# LO SOLUTION             -4.2757D+00   $ (NH=400)
# LO SOLUTION             -4.2743D+00   $ (NH=800)
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
        pb.pbclass = "C-LOR2-AN-V-V"
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

    elseif action == "eSQR"

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

