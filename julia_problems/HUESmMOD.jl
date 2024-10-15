function HUESmMOD(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HUESmMOD
#    *********
# 
#    Source: An inverse problem from astronomy,
#    reformulated as a convex quadratic program by
#    S. P. Hestis, SIAM Review 34 (1992) pp. 642-647.
# 
#    SIF input: Nick Gould, January 1993.
#    improvements by: Ruediger Franke (Ruediger.Franke@RZ.TU-Ilmenau.DE)
# 
#    classification = "C-QLR2-MN-V-V"
# 
#    Number of variables
# 
#       Alternative values for the SIF file parameters:
# IE K                   10             $-PARAMETER
# IE K                   100            $-PARAMETER
# IE K                   1000           $-PARAMETER    original value
# IE K                   5000           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HUESmMOD"

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
            v_["K"] = Int64(10);  #  SIF file default value
        else
            v_["K"] = Int64(args[1]);
        end
# IE K                   10000          $-PARAMETER
        v_["1"] = 1
        v_["RANGE"] = 1.0
        v_["3.0"] = 3.0
        v_["5.0"] = 5.0
        v_["RK"] = Float64(v_["K"])
        v_["DELTAX"] = v_["RANGE"]/v_["RK"]
        v_["DELTAX2"] = v_["DELTAX"]*v_["DELTAX"]
        v_["DELTAX3"] = v_["DELTAX2"]*v_["DELTAX"]
        v_["DELTAX5"] = v_["DELTAX3"]*v_["DELTAX2"]
        v_["DELTAX3/3"] = v_["DELTAX3"]/v_["3.0"]
        v_["DELTAX5/5"] = v_["DELTAX5"]/v_["5.0"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["K"])
            iv,ix_,_ = s2mpj_ii("M"*string(I),ix_)
            arrset(pb.xnames,iv,"M"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["K"])
            ig,ig_,_ = s2mpj_ii("OBJ"*string(I),ig_)
            arrset(gtype,ig,"<>")
        end
        for I = Int64(v_["1"]):Int64(v_["K"])
            v_["I-1"] = -1+I
            v_["RI"] = Float64(I)
            v_["RI2"] = v_["RI"]*v_["RI"]
            v_["RI3"] = v_["RI2"]*v_["RI"]
            v_["RI5"] = v_["RI3"]*v_["RI2"]
            v_["RI-1"] = Float64(v_["I-1"])
            v_["RI-12"] = v_["RI-1"]*v_["RI-1"]
            v_["RI-13"] = v_["RI-12"]*v_["RI-1"]
            v_["RI-15"] = v_["RI-13"]*v_["RI-12"]
            v_["DIFF3"] = v_["RI3"]-v_["RI-13"]
            v_["DIFF5"] = v_["RI5"]-v_["RI-15"]
            v_["COEFF1"] = v_["DIFF3"]*v_["DELTAX3/3"]
            v_["COEFF2"] = v_["DIFF5"]*v_["DELTAX5/5"]
            ig,ig_,_ = s2mpj_ii("E1",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E1")
            iv = ix_["M"*string(I)]
            pbm.A[ig,iv] += Float64(v_["COEFF1"])
            ig,ig_,_ = s2mpj_ii("E2",ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"E2")
            iv = ix_["M"*string(I)]
            pbm.A[ig,iv] += Float64(v_["COEFF2"])
        end
        v_["RK"] = Float64(v_["K"])
        v_["1/RK"] = 1.0/v_["RK"]
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
        pbm.gconst[ig_["E1"]] = Float64(1835.2)
        pbm.gconst[ig_["E2"]] = Float64(909.8)
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        pb.y0 = fill(Float64(1.0),pb.m)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"U1")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["K"])
            ename = "E"*string(I)
            ie,ie_,newelt = s2mpj_ii(ename,ie_)
            if newelt > 0
                arrset(pbm.elftype,ie,"eSQ")
                arrset(ielftype,ie,iet_["eSQ"])
            end
            vname = "M"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="U1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["K"])
            ig = ig_["OBJ"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["1/RK"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0
#    Solution
# LO SOLTN               0.0
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
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
        pb.pbclass = "C-QLR2-MN-V-V"
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

