function CLNLBEAM(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CLNLBEAM
#    *********
# 
#    An optimal control version of the CLamped NonLinear BEAM problem.
#    The energy of a beam of length 1 compressed by a force P is to be
#    minimized.  The control variable is the derivative of the deflection angle.
# 
#    The problem is discretized using the trapezoidal rule. It is non-convex.
# 
#    Source:
#    H. Maurer and H.D. Mittelman,
#    "The non-linear beam via optimal control with bound state variables",
#    Optimal Control Applications and Methods 12, pp. 19-31, 1991.
# 
#    SIF input: Ph. Toint, Nov 1993.
# 
#    classification = "C-OOR2-MN-V-V"
# 
#    Discretization: specify the number of interior points + 1
# 
#       Alternative values for the SIF file parameters:
# IE NI                  10             $-PARAMETER n=33, m=20
# IE NI                  50             $-PARAMETER n=153, m=100
# IE NI                  100            $-PARAMETER n=303, m=200
# IE NI                  500            $-PARAMETER n=1503, m=1000
# IE NI                  1000           $-PARAMETER n=3003, m=2000 original value
# IE NI                  2000           $-PARAMETER n=6003, m=4000
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CLNLBEAM"

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
            v_["NI"] = Int64(10);  #  SIF file default value
        else
            v_["NI"] = Int64(args[1]);
        end
# IE NI                  5000           $-PARAMETER n=15003, m=10000
        if nargin<2
            v_["ALPHA"] = Float64(350.0);  #  SIF file default value
        else
            v_["ALPHA"] = Float64(args[2]);
        end
        v_["RNI"] = Float64(v_["NI"])
        v_["NI-1"] = -1+v_["NI"]
        v_["H"] = 1.0/v_["RNI"]
        v_["H/4"] = 0.25*v_["H"]
        v_["H/2"] = 0.5*v_["H"]
        v_["AH"] = v_["ALPHA"]*v_["H"]
        v_["AH/2"] = 0.5*v_["AH"]
        v_["-H/2"] = -0.5*v_["H"]
        v_["0"] = 0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NI"])
            iv,ix_,_ = s2mpj_ii("T"*string(I),ix_)
            arrset(pb.xnames,iv,"T"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["NI"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        for I = Int64(v_["0"]):Int64(v_["NI"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("ENERGY",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            v_["I+1"] = 1+I
            ig,ig_,_ = s2mpj_ii("EX"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"EX"*string(I))
            iv = ix_["X"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ET"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ET"*string(I))
            iv = ix_["T"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(1.0)
            iv = ix_["T"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            iv = ix_["U"*string(Int64(v_["I+1"]))]
            pbm.A[ig,iv] += Float64(v_["-H/2"])
            iv = ix_["U"*string(I)]
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
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        for I = Int64(v_["0"]):Int64(v_["NI"])
            pb.xlower[ix_["X"*string(I)]] = -0.05
            pb.xupper[ix_["X"*string(I)]] = 0.05
        end
        for I = Int64(v_["0"]):Int64(v_["NI"])
            pb.xlower[ix_["T"*string(I)]] = -1.0
            pb.xupper[ix_["T"*string(I)]] = 1.0
        end
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["NI"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["NI"]))]] = 0.0
        pb.xlower[ix_["T"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["T"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["T"*string(Int64(v_["NI"]))]] = 0.0
        pb.xupper[ix_["T"*string(Int64(v_["NI"]))]] = 0.0
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        for I = Int64(v_["0"]):Int64(v_["NI"])
            v_["RI"] = Float64(I)
            v_["TT"] = v_["RI"]*v_["H"]
            v_["CTT"] = cos(v_["TT"])
            v_["SCTT"] = 0.05*v_["CTT"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["SCTT"])
            pb.x0[ix_["T"*string(I)]] = Float64(v_["SCTT"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eCOS", iet_)
        loaset(elftv,it,1,"T")
        it,iet_,_ = s2mpj_ii( "eSIN", iet_)
        loaset(elftv,it,1,"T")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"U")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["0"]):Int64(v_["NI"])
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eCOS")
            arrset(ielftype,ie,iet_["eCOS"])
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "S"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSIN")
            arrset(ielftype,ie,iet_["eSIN"])
            vname = "T"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="T",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "USQ"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="U",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["0"]):Int64(v_["NI-1"])
            v_["I+1"] = 1+I
            ig = ig_["EX"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["S"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["-H/2"]))
            ig = ig_["ENERGY"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["USQ"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["USQ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(Int64(v_["I+1"]))])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["AH/2"]))
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["C"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["AH/2"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           345.0301196587
# LO SOLTN(50)           344.8673691861
# LO SOLTN(100)          344.8801831150
# LO SOLTN(500)          344.8748539754
# LO SOLTN(1000)         344.8788169123
# LO SOLTN(5000)         
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

    elseif action == "eCOS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CC = cos(EV_[1])
        SS = sin(EV_[1])
        f_   = CC
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -SS
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -CC
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eSIN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        CC = cos(EV_[1])
        SS = sin(EV_[1])
        f_   = SS
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = CC
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = -SS
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

