function OPTMASS(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : OPTMASS
#    *********
# 
#    A constrained optimal control problem
#    adapted from Gawande and Dunn
# 
#    The problem is that of a particle of unit mass moving on a
#    frictionless plane under the action of a controlling force whose
#    magnitude may not exceed unity. At time=0, the particle moves through
#    the origin of the plane in the direction of the positive x-axis with
#    speed SPEED.  The cost function incorporates two conflicting control
#    objectives, namely: maximization of the particle's final (at time=1)
#    distance from the origin and minimization of its final speed.  By
#    increasing the  value of the penalty constant PEN, more stress can be
#    placed on the latter objective.
# 
#    Gawande and Dunn originally use a starting point (in the control
#    only) that is much closer to the solution than the one chosen
#    here.
# 
#    Source:
#    M. Gawande and J. Dunn,
#    "A Projected Newton Method in a Cartesian Product of Balls",
#    JOTA 59(1): 59-69, 1988.
# 
#    SIF input: Ph. Toint, June 1990.
# 
#    classification = "C-QQR2-AN-V-V"
# 
#    Number of discretization steps in the time interval
#    The number of variables is 6 * (N + 2) -2 , 4 of which are fixed.
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER n = 70    original value
# IE N                   100            $-PARAMETER n = 610
# IE N                   200            $-PARAMETER n = 1210
# IE N                   500            $-PARAMETER n = 3010
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "OPTMASS"

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
# IE N                   1000           $-PARAMETER n = 6010
# IE N                   5000           $-PARAMETER n = 30010
        v_["SPEED"] = 0.01
        v_["PEN"] = 0.335
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        v_["N+1"] = 1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["1/N"] = 1.0/v_["RN"]
        v_["-1/N"] = -1.0*v_["1/N"]
        v_["1/N2"] = v_["1/N"]*v_["1/N"]
        v_["-1/2N2"] = -0.5*v_["1/N2"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["2"])
                iv,ix_,_ = s2mpj_ii("X"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"X"*string(J)*","*string(I))
                iv,ix_,_ = s2mpj_ii("V"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"V"*string(J)*","*string(I))
                iv,ix_,_ = s2mpj_ii("F"*string(J)*","*string(I),ix_)
                arrset(pb.xnames,iv,"F"*string(J)*","*string(I))
            end
        end
        for J = Int64(v_["1"]):Int64(v_["2"])
            iv,ix_,_ = s2mpj_ii("X"*string(J)*","*string(Int64(v_["N+1"])),ix_)
            arrset(pb.xnames,iv,"X"*string(J)*","*string(Int64(v_["N+1"])))
            iv,ix_,_ = s2mpj_ii("V"*string(J)*","*string(Int64(v_["N+1"])),ix_)
            arrset(pb.xnames,iv,"V"*string(J)*","*string(Int64(v_["N+1"])))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("F",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            v_["I-1"] = -1+I
            for J = Int64(v_["1"]):Int64(v_["2"])
                ig,ig_,_ = s2mpj_ii("A"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"A"*string(J)*","*string(I))
                iv = ix_["X"*string(J)*","*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["X"*string(J)*","*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["V"*string(J)*","*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(v_["-1/N"])
                iv = ix_["F"*string(J)*","*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(v_["-1/2N2"])
                ig,ig_,_ = s2mpj_ii("B"*string(J)*","*string(I),ig_)
                arrset(gtype,ig,"==")
                arrset(pb.cnames,ig,"B"*string(J)*","*string(I))
                iv = ix_["V"*string(J)*","*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["V"*string(J)*","*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(-1.0)
                iv = ix_["F"*string(J)*","*string(Int64(v_["I-1"]))]
                pbm.A[ig,iv] += Float64(v_["-1/N"])
            end
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("C"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"C"*string(I))
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
        for I = Int64(v_["0"]):Int64(v_["N"])
            pbm.gconst[ig_["C"*string(I)]] = Float64(1.0)
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["V"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
              v_["SPEED"])
        pb.xupper[ix_["V"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
              v_["SPEED"])
        pb.xlower[ix_["V"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["V"*string(Int64(v_["2"]))*","*string(Int64(v_["0"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.0),pb.n)
        if haskey(ix_,"V"*string(Int64(v_["1"]))*","*string(Int64(v_["0"])))
            pb.x0[ix_["V"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))]]  = (
                  Float64(v_["SPEED"]))
        else
            pb.y0[findfirst(x->x==ig_["V"*string(Int64(v_["1"]))*","*string(Int64(v_["0"]))],pbm.congrps)] = Float64(v_["SPEED"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        ename = "O1"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X"*string(Int64(v_["1"]))*","*string(Int64(v_["N+1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "O2"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "X"*string(Int64(v_["2"]))*","*string(Int64(v_["N+1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "O3"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "V"*string(Int64(v_["1"]))*","*string(Int64(v_["N+1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        ename = "O4"
        ie,ie_,_  = s2mpj_ii(ename,ie_)
        arrset(pbm.elftype,ie,"eSQ")
        arrset(ielftype,ie,iet_["eSQ"])
        vname = "V"*string(Int64(v_["2"]))*","*string(Int64(v_["N+1"]))
        iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
        posev = findfirst(x->x=="X",elftv[ielftype[ie]])
        loaset(pbm.elvar,ie,posev,iv)
        for I = Int64(v_["0"]):Int64(v_["N"])
            for J = Int64(v_["1"]):Int64(v_["2"])
                ename = "D"*string(J)*","*string(I)
                ie,ie_,_  = s2mpj_ii(ename,ie_)
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
                vname = "F"*string(J)*","*string(I)
                iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.0))
                posev = findfirst(x->x=="X",elftv[ielftype[ie]])
                loaset(pbm.elvar,ie,posev,iv)
            end
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["F"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["O1"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["O2"])
        loaset(pbm.grelw,ig,posel,Float64(-1.0))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["O3"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["PEN"]))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["O4"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(v_["PEN"]))
        for I = Int64(v_["0"]):Int64(v_["N"])
            ig = ig_["C"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["1"]))*","*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(Int64(v_["2"]))*","*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           -0.04647
# LO SOLTN(100)          ???
# LO SOLTN(200)          ???
# LO SOLTN(500)          ???
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-QQR2-AN-V-V"
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

