function PENLT2NE(action::String,args::Union{Any}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PENLT2NE
#    --------
# 
#    The second penalty function
# 
#    This is a nonlinear least-squares problem with M=2*N groups.
#     Group 1 is linear.
#     Groups 2 to N use 2 nonlinear elements.
#     Groups N+1 to M-1 use 1 nonlinear element.
#     Group M uses N nonlinear elements.
#    The Hessian matrix is dense. This is a nonlinear equation version
#    of problem PENALTY2.
# 
#    Source:  Problem 24 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
# 
#    See also Buckley#112 (p. 80)
# 
#    SIF input: Ph. Toint, Dec 1989.
#    Modification as a set of nonlinear equations: Nick Gould, Oct 2015.
# 
#    classification = "C-CNOR2-AN-V-V"
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE N                   4              $-PARAMETER     original value
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 21 VI 2025
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PENLT2NE"
    if ( !isdefined(@__MODULE__, :s2mpj_ii) )
        error( "Please include(\"s2mpjlib.jl\") using \"s2mpjlib.jl\" from the S2MPJ distribution before calling PENLT2NE.")
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
# IE N                   10             $-PARAMETER
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   200            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
        v_["A"] = 0.00001
        v_["B"] = 1.0
        v_["1"] = 1
        v_["2"] = 2
        v_["N+1"] = 1+v_["N"]
        v_["M"] = v_["N"]+v_["N"]
        v_["M-1"] = -1+v_["M"]
        v_["EM1/10"] = exp(-0.1)
        v_["1/A"] = 1.0/v_["A"]
        v_["1/B"] = 1.0/v_["B"]
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
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G"*string(Int64(v_["1"])))
        push!(irA,ig)
        push!(icA,ix_["X"*string(Int64(v_["1"]))])
        push!(valA,Float64(1.0))
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["1"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G"*string(Int64(v_["1"])))
        arrset(pbm.gscale,ig,Float64(v_["1/B"]))
        for I = Int64(v_["2"]):Int64(v_["M-1"])
            ig,ig_,_ = s2mpj_ii("G"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"G"*string(I))
            arrset(pbm.gscale,ig,Float64(v_["1/A"]))
        end
        ig,ig_,_ = s2mpj_ii("G"*string(Int64(v_["M"])),ig_)
        arrset(gtype,ig,"==")
        arrset(pb.cnames,ig,"G"*string(Int64(v_["M"])))
        arrset(pbm.gscale,ig,Float64(v_["1/B"]))
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
        pbm.gconst[ig_["G1"]] = Float64(0.2)
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["RI-1"] = Float64(v_["I-1"])
            v_["I/10"] = 0.1*v_["RI"]
            v_["I-1/10"] = 0.1*v_["RI-1"]
            v_["EI/10"] = exp(v_["I/10"])
            v_["EI-1/10"] = exp(v_["I-1/10"])
            v_["YI"] = v_["EI/10"]+v_["EI-1/10"]
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["YI"])
        end
        for I = Int64(v_["N+1"]):Int64(v_["M-1"])
            pbm.gconst[ig_["G"*string(I)]] = Float64(v_["EM1/10"])
        end
        pbm.gconst[ig_["G"*string(Int64(v_["M"]))]] = Float64(1.0)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.5),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE10", iet_)
        loaset(elftv,it,1,"V")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"V")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["2"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE10")
            arrset(ielftype,ie,iet_["eE10"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE10")
            arrset(ielftype,ie,iet_["eE10"])
            vname = "X"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["N+1"]):Int64(v_["M-1"])
            v_["-N"] = -1*v_["N"]
            v_["I-N"] = I+v_["-N"]
            v_["I-N+1"] = 1+v_["I-N"]
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE10")
            arrset(ielftype,ie,iet_["eE10"])
            vname = "X"*string(Int64(v_["I-N+1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            ename = "D"*string(J)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "X"*string(J)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.5))
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["2"]):Int64(v_["N"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        for I = Int64(v_["N+1"]):Int64(v_["M-1"])
            ig = ig_["G"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        for J = Int64(v_["1"]):Int64(v_["N"])
            v_["-J"] = -1*J
            v_["N-J"] = v_["N"]+v_["-J"]
            v_["N-J+1"] = 1+v_["N-J"]
            v_["WI"] = Float64(v_["N-J+1"])
            ig = ig_["G"*string(Int64(v_["M"]))]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["D"*string(J)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["WI"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN(4)            9.37629D-6
# LO SOLTN(10)           2.93660D-4
# LO SOLTN(50)           4.29609813
# LO SOLTN(100)          97096.0840
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
        pb.pbclass = "C-CNOR2-AN-V-V"
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

    elseif action == "eE10"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        EXPA = exp(0.1*EV_[1])
        f_   = EXPA
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 0.1*EXPA
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 0.01*EXPA
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

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

