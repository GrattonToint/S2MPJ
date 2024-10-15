function TRAINF(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : TRAINF
#    *********
# 
#    This is an optimal control problem.
#    The problem is to minimize the energy spent to move a train 
#    from the beginning of a flat track to its end in a given time.  The train
#    is slowed down by some drag (assumed to be quadratic in the the velocity).
#    The control variables are the acceleration force (UA) and the braking
#    force (UB) applied on the train.
# 
#    Source:
#    J. Kautsky and N. K. Nichols,
#    "OTEP-2: Optimal Train Energy Programme, mark 2",
#    Numerical Analysis Report NA/4/83,
#    Department of Mathematics, University of Reading, 1983.
# 
#    SIF input: N. Nichols and Ph. Toint, April 1993
# 
#    classification = "C-QQR2-MN-V-V"
# 
#    Problem variants
# 
#       Alternative values for the SIF file parameters:
# RE TIME                4.8            $-PARAMETER  travel time
# RE LENGTH              6.0            $-PARAMETER  length of track
# 
# RE TIME                2.0            $-PARAMETER  travel time
# RE LENGTH              2.0            $-PARAMETER  length of track
# 
# RE TIME                1.5            $-PARAMETER  travel time
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "TRAINF"

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
            v_["TIME"] = Float64(1.5);  #  SIF file default value
        else
            v_["TIME"] = Float64(args[1]);
        end
# RE LENGTH              2.0            $-PARAMETER  length of track
        if nargin<2
            v_["LENGTH"] = Float64(2);  #  SIF file default value
        else
            v_["LENGTH"] = Float64(args[2]);
        end
# IE N                   11             $-PARAMETER
# IE N                   51             $-PARAMETER
# IE N                   101            $-PARAMETER     original value
# IE N                   201            $-PARAMETER
# IE N                   501            $-PARAMETER
# IE N                   1001           $-PARAMETER
        if nargin<3
            v_["N"] = Int64(11);  #  SIF file default value
        else
            v_["N"] = Int64(args[3]);
        end
# IE N                   5001           $-PARAMETER
# IE N                   10001          $-PARAMETER
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = v_["TIME"]/v_["RN"]
        v_["H/2"] = 0.5*v_["H"]
        v_["-H"] = -1.0*v_["H"]
        v_["-H/2"] = -1.0*v_["H/2"]
        v_["UAMAX"] = 10.0
        v_["UBMIN"] = -2.0
        v_["VMAX"] = 10.0
        v_["A"] = 0.3
        v_["B"] = 0.14
        v_["C"] = 0.16
        v_["0"] = 0
        v_["1"] = 1
        v_["BH/2"] = v_["B"]*v_["H/2"]
        v_["1+BH/2"] = 1.0+v_["BH/2"]
        v_["BH/2-1"] = -1.0+v_["BH/2"]
        v_["-AH"] = v_["A"]*v_["-H"]
        v_["LENGTH/N"] = v_["LENGTH"]/v_["RN"]
        v_["CH/2"] = v_["C"]*v_["H/2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("V"*string(I),ix_)
            arrset(pb.xnames,iv,"V"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("UA"*string(I),ix_)
            arrset(pb.xnames,iv,"UA"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("UB"*string(I),ix_)
            arrset(pb.xnames,iv,"UB"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("ENERGY",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("XEQ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"XEQ"*string(I))
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["V"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["V"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            ig,ig_,_ = s2mpj_ii("VEQ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"VEQ"*string(I))
            iv = ix_["V"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["1+BH/2"])
            iv = ix_["V"*string(I)]
            pbm.A[ig,iv] += Float64(v_["BH/2-1"])
            iv = ix_["UA"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UA"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UB"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["UB"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
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
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            pbm.gconst[ig_["VEQ"*string(I)]] = Float64(v_["-AH"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["V"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["UA"*string(Int64(v_["0"]))]] = v_["UAMAX"]
        pb.xupper[ix_["UA"*string(Int64(v_["0"]))]] = v_["UAMAX"]
        pb.xlower[ix_["UB"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["UB"*string(Int64(v_["0"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            pb.xlower[ix_["X"*string(I)]] = -Inf
            pb.xupper[ix_["X"*string(I)]] = +Inf
            pb.xlower[ix_["V"*string(I)]] = -Inf
            pb.xupper[ix_["V"*string(I)]] = +Inf
            pb.xlower[ix_["UA"*string(I)]] = 0.0
            pb.xupper[ix_["UA"*string(I)]] = v_["UAMAX"]
            pb.xlower[ix_["UB"*string(I)]] = v_["UBMIN"]
            pb.xupper[ix_["UB"*string(I)]] = 0.0
        end
        pb.xlower[ix_["X"*string(Int64(v_["N"]))]] = v_["LENGTH"]
        pb.xupper[ix_["X"*string(Int64(v_["N"]))]] = v_["LENGTH"]
        pb.xlower[ix_["V"*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["UA"*string(Int64(v_["N"]))]] = 0.0
        pb.xupper[ix_["UA"*string(Int64(v_["N"]))]] = 0.0
        pb.xlower[ix_["UB"*string(Int64(v_["N"]))]] = v_["UBMIN"]
        pb.xupper[ix_["UB"*string(Int64(v_["N"]))]] = v_["UBMIN"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["V"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["UA"*string(Int64(v_["0"]))]] = Float64(v_["UAMAX"])
        pb.x0[ix_["UB"*string(Int64(v_["0"]))]] = Float64(0.0)
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["RI"] = Float64(I)
            v_["PI"] = v_["LENGTH/N"]*v_["RI"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["PI"])
            pb.x0[ix_["V"*string(I)]] = Float64(v_["LENGTH/N"])
            pb.x0[ix_["UA"*string(I)]] = Float64(0.0)
            pb.x0[ix_["UB"*string(I)]] = Float64(0.0)
        end
        pb.x0[ix_["X"*string(Int64(v_["N"]))]] = Float64(v_["LENGTH"])
        pb.x0[ix_["V"*string(Int64(v_["N"]))]] = Float64(0.0)
        pb.x0[ix_["UA"*string(Int64(v_["N"]))]] = Float64(0.0)
        pb.x0[ix_["UB"*string(Int64(v_["N"]))]] = Float64(v_["UBMIN"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"UU")
        loaset(elftv,it,2,"VV")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"VVV")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["N"])
            v_["I+1"] = 1+I
            ename = "VISQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "V"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VVV",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ename = "UV"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "UA"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="UU",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "V"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="VV",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig = ig_["VEQ"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["VISQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["CH/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["VISQ"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["CH/2"]))
        end
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            ig = ig_["ENERGY"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["UV"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
# LO SOLUTION            3.09751881012
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
        pb.pbclass = "C-QQR2-MN-V-V"
        pbm.objderlvl = 2
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2]
        pb.conderlvl  = pbm.conderlvl;
        return pb, pbm

# ********************
#  SET UP THE GROUPS *
#  ROUTINE           *
# ********************

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

