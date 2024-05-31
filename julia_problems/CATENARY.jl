function CATENARY(action,args...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    An erroneous but interesting version of the classical problem of the hanging
#    catenary.  Here the catenary consists of N+1 beams of length BL,
#    with the first beam fixed at the origin and the final bean fixed
#    at a fraction FRACT of the total length of all beams.
# 
#    The correct version of the same problem is given by problem CATENA.SIF
# 
#    The problem is non-convex.
# 
#    Source: 
#    K. Veselic,
#    "De forma catenarum in campo gravitatis pendentium",
#    Klasicna Gimnazija u Zagrebu, Zagreb, 1987.
# 
#    SIF input: Ph. L. Toint, May 1993.
# 
#    classification = "LQR2-AY-V-V"
# 
#    Number of beams = N+1 ; the number of variables is 3*(N+2)
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CATENARY"

    if action == "setup"
        pbm          = PBM(name)
        pb           = PB(name)
        pb.sifpbname = "CATENARY"
        nargin       = length(args)
        pbm.call     = eval( Meta.parse( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        if nargin<1
            v_["N+1"] = Int64(4);  #  SIF file default value
        else
            v_["N+1"] = Int64(args[1]);
        end
#       Alternative values for the SIF file parameters:
# IE N+1                 10             $-PARAMETER n = 33
# IE N+1                 32             $-PARAMETER n = 99
# IE N+1                 166            $-PARAMETER n = 501     original value
# IE N+1                 1000           $-PARAMETER n = 3003
        v_["GAMMA"] = 9.81
        v_["TMASS"] = 500.0
        v_["BL"] = 1.0
        v_["FRACT"] = 0.6
        v_["0"] = 0
        v_["1"] = 1
        v_["N"] = -1+v_["N+1"]
        v_["RN+1"] = Float64(v_["N+1"])
        v_["TMP"] = v_["BL"]*v_["RN+1"]
        v_["LENGTH"] = v_["TMP"]*v_["FRACT"]
        v_["MASS"] = v_["TMASS"]/v_["RN+1"]
        v_["MG"] = v_["MASS"]*v_["GAMMA"]
        v_["MG/2"] = 0.5*v_["MG"]
        v_["BLSQ"] = v_["BL"]*v_["BL"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        xscale  = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N+1"])
            iv,ix_,_ = s2x_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2x_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2x_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Y"*string(Int64(v_["0"]))]
        pbm.A[ig,iv] += Float64(v_["MG/2"])
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2x_ii("OBJ",ig_)
            arrset(gtype,ig,"<>")
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += Float64(v_["MG"])
        end
        ig,ig_,_ = s2x_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        iv = ix_["Y"*string(Int64(v_["N+1"]))]
        pbm.A[ig,iv] += Float64(v_["MG/2"])
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            ig,ig_,_ = s2x_ii("C"*string(I),ig_)
            arrset(gtype,ig,"==")
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
        pbm.congrps = findall(x->x!="<>",gtype)
        pb.nob = ngrp-pb.m
        pbm.objgrps = findall(x->x=="<>",gtype)
        #%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(Float64,ngrp)
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            pbm.gconst[ig_["C"*string(I)]] = Float64(v_["BLSQ"])
        end
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X0"]] = 0.0
        pb.xupper[ix_["X0"]] = 0.0
        pb.xlower[ix_["Y0"]] = 0.0
        pb.xupper[ix_["Y0"]] = 0.0
        pb.xlower[ix_["Z0"]] = 0.0
        pb.xupper[ix_["Z0"]] = 0.0
        pb.xlower[ix_["X"*string(Int64(v_["N+1"]))]] = v_["LENGTH"]
        pb.xupper[ix_["X"*string(Int64(v_["N+1"]))]] = v_["LENGTH"]
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        v_["TMP"] = v_["LENGTH"]/v_["RN+1"]
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            v_["RI"] = Float64(I)
            v_["VAL"] = v_["TMP"]*v_["RI"]
            pb.x0[ix_["X"*string(I)]] = Float64(v_["TMP"])
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2x_ii( "eISQ", iet_)
        loaset(elftv,it,1,"V")
        loaset(elftv,it,2,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            v_["I-1"] = -1+I
            ename = "EX"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eISQ")
            arrset(ielftype, ie, iet_["eISQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EY"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eISQ")
            arrset(ielftype, ie, iet_["eISQ"])
            vname = "Y"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EZ"*string(I)
            ie,ie_,_  = s2x_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eISQ")
            arrset(ielftype, ie, iet_["eISQ"])
            vname = "X"*string(I)
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2x_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="W",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N+1"])
            ig = ig_["C"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EX"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EY"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EZ"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "LQR2-AY-V-V"
        return pb, pbm

    #%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    elseif action == "eISQ"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        U_ = zeros(Float64,1,2)
        IV_ =  zeros(Float64,1)
        U_[1,1] = U_[1,1]+1
        U_[1,2] = U_[1,2]-1
        IV_[1] = dot(U_[1,:],EV_)
        f_   = IV_[1]*IV_[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = IV_[1]+IV_[1]
            g_ =  U_'*g_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_[1,1] = 2.0
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

    #%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    elseif action in  ["fx","fgx","fgHx","cx","cJx","cJHx","cIx","cIJx","cIJHx","cIJxv","fHxv","cJxv","Lxy","Lgxy","LgHxy","LIxy","LIgxy","LIgHxy","LHxyv","LIHxyv"]

        pbm = args[1]
        if pbm.name == name
            pbm.has_globs = [0,0]
            return s2x_eval(action,args...)
        else
            println("ERROR: please run "*name*" with action = setup")
            return ntuple(i->undef,args[end])
        end

    else
        println("ERROR: unknown action "*action*" requested from "*name*"%s.jl")
        return ntuple(i->undef,args[end])
    end

end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
