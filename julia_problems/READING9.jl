function READING9(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : READING9
#    *********
# 
#    A nonlinear optimal control problem from Nancy Nichols
#    with a given initial condition.
#    This problem arises in tide modelling.
# 
#    Source: a variant upon a problem in
#    S. Lyle and N.K. Nichols,
#    "Numerical Methods for Optimal Control Problems with State Constraints",
#    Numerical Analysis Report 8/91, Dept of Mathematics, 
#    University of Reading, UK.
# 
#    SIF input: Nick Gould and Ph. Toint, March 1995
# 
#    classification = "C-OOR2-MN-V-V"
# 
#    Number of discretized points in [0,1] - 1
# 
#       Alternative values for the SIF file parameters:
# IE N+1                 3              $-PARAMETER n=6, m=2
# IE N+1                 51             $-PARAMETER n=102, m=50
# IE N+1                 101            $-PARAMETER n=202, m=100
# IE N+1                 201            $-PARAMETER n=402, m=200
# IE N+1                 501            $-PARAMETER n=1002, m=500  original value
# IE N+1                 1001           $-PARAMETER n=2002, m=1000
# IE N+1                 5001           $-PARAMETER n=10002, m= 5000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "READING9"

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
            v_["N+1"] = Int64(3);  #  SIF file default value
        else
            v_["N+1"] = Int64(args[1]);
        end
        v_["N"] = -1+v_["N+1"]
        v_["N-1"] = -1+v_["N"]
        v_["RN"] = Float64(v_["N"])
        v_["H"] = 1.0/v_["RN"]
        v_["-H/2"] = -0.5*v_["H"]
        v_["K1"] = 0.07716
        v_["K2"] = 0.07716
        v_["K1H"] = v_["K1"]*v_["H"]
        v_["K1H+1"] = 1.0+v_["K1H"]
        v_["-K1H-1"] = -1.0*v_["K1H+1"]
        v_["K2H"] = v_["K2"]*v_["H"]
        v_["0"] = 0
        v_["1"] = 1
        v_["2"] = 2
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("P"*string(I),ix_)
            arrset(pb.xnames,iv,"P"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("S"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"S"*string(I))
            iv = ix_["P"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["P"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-K1H-1"])
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
            v_["RI"] = Float64(I)
            v_["T"] = v_["RI"]*v_["H"]
            v_["SINT"] = sin(v_["T"])
            v_["HSINT"] = v_["H"]*v_["SINT"]
            pbm.gconst[ig_["S"*string(I)]] = Float64(v_["HSINT"])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["0"]):Int64(v_["N"])
            pb.xupper[ix_["U"*string(I)]] = 1.0
            pb.xlower[ix_["P"*string(I)]] = -Inf
            pb.xupper[ix_["P"*string(I)]] = +Inf
        end
        pb.xlower[ix_["P"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["P"*string(Int64(v_["0"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(0.2),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "ePROD", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"P")
        it,iet_,_ = s2mpj_ii( "ePROD2", iet_)
        loaset(elftv,it,1,"U")
        loaset(elftv,it,2,"P")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["N"])
            ename = "OE"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD2")
            arrset(ielftype,ie,iet_["ePROD2"])
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.2))
            posev = findfirst(x->x=="P",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.2))
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            ename = "CE"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"ePROD")
            arrset(ielftype,ie,iet_["ePROD"])
            vname = "P"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.2))
            posev = findfirst(x->x=="P",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(0.2))
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N-1"])
            v_["I+1"] = 1+I
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["OE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["OE"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            ig = ig_["S"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["CE"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["K2H"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# XL SOLUTION            -4.41677D-02   $ (n=500)
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
        pb.pbclass = "C-OOR2-MN-V-V"
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

    elseif action == "ePROD2"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_   = (EV_[1]*EV_[2])^2
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*EV_[1]*EV_[2]^2
            g_[2] = 2.0*EV_[2]*EV_[1]^2
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*EV_[2]^2
                H_[1,2] = 4.0*EV_[1]*EV_[2]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[1]^2
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

