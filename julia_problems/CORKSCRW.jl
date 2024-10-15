function CORKSCRW(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CORKSCRW
#    *********
# 
#    A nonlinear optimal control problem with both state- and
#    control constraints.
#    The problem is to control (using an applied force of limited
#    magnitude) a mass moving in the 3D space, such that its
#    trajectory lies within a prescribed distance TOL of the
#    corkscreww-like curve defined by
#               y = sin(x), z = cos(x),
#    and such that it stops at a given abscissa XT in minimum time.
#    The mass is initially stationary at (0,0,1).
# 
#    Source:
#    Ph. Toint, private communication.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-SOR2-AN-V-V"
# 
#    Number of time intervals
#    The number of variables is 9T+6, of which 9 are fixed.
# 
#       Alternative values for the SIF file parameters:
# IE T                   10             $-PARAMETER n = 96     original value
# IE T                   50             $-PARAMETER n = 456
# IE T                   100            $-PARAMETER n = 906
# IE T                   500            $-PARAMETER n = 4506
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CORKSCRW"

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
            v_["T"] = Int64(10);  #  SIF file default value
        else
            v_["T"] = Int64(args[1]);
        end
# IE T                   1000           $-PARAMETER n = 9006
        if nargin<2
            v_["XT"] = Float64(10.0);  #  SIF file default value
        else
            v_["XT"] = Float64(args[2]);
        end
        if nargin<3
            v_["MASS"] = Float64(0.37);  #  SIF file default value
        else
            v_["MASS"] = Float64(args[3]);
        end
        if nargin<4
            v_["TOL"] = Float64(0.1);  #  SIF file default value
        else
            v_["TOL"] = Float64(args[4]);
        end
        v_["0"] = 0
        v_["1"] = 1
        v_["RT"] = Float64(v_["T"])
        v_["T+1"] = 1.0+v_["RT"]
        v_["H"] = v_["XT"]/v_["RT"]
        v_["1/H"] = 1.0/v_["H"]
        v_["-1/H"] = -1.0*v_["1/H"]
        v_["M/H"] = v_["MASS"]/v_["H"]
        v_["-M/H"] = -1.0*v_["M/H"]
        v_["TOLSQ"] = v_["TOL"]*v_["TOL"]
        v_["XTT+1"] = v_["XT"]*v_["T+1"]
        v_["W"] = 0.5*v_["XTT+1"]
        for I = Int64(v_["1"]):Int64(v_["T"])
            v_["RI"] = Float64(I)
            v_["TI"] = v_["RI"]*v_["H"]
            v_["W/T"*string(I)] = v_["W"]/v_["TI"]
        end
        v_["FMAX"] = v_["XT"]/v_["RT"]
        v_["-FMAX"] = -1.0*v_["FMAX"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["T"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
            iv,ix_,_ = s2mpj_ii("Z"*string(I),ix_)
            arrset(pb.xnames,iv,"Z"*string(I))
            iv,ix_,_ = s2mpj_ii("VX"*string(I),ix_)
            arrset(pb.xnames,iv,"VX"*string(I))
            iv,ix_,_ = s2mpj_ii("VY"*string(I),ix_)
            arrset(pb.xnames,iv,"VY"*string(I))
            iv,ix_,_ = s2mpj_ii("VZ"*string(I),ix_)
            arrset(pb.xnames,iv,"VZ"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["T"])
            iv,ix_,_ = s2mpj_ii("UX"*string(I),ix_)
            arrset(pb.xnames,iv,"UX"*string(I))
            iv,ix_,_ = s2mpj_ii("UY"*string(I),ix_)
            arrset(pb.xnames,iv,"UY"*string(I))
            iv,ix_,_ = s2mpj_ii("UZ"*string(I),ix_)
            arrset(pb.xnames,iv,"UZ"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["T"])
            ig,ig_,_ = s2mpj_ii("OX"*string(I),ig_)
            arrset(gtype,ig,"<>")
            arrset(pbm.gscale,ig,Float64(v_["W/T"*string(I)]))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(1.0)
        end
        for I = Int64(v_["1"]):Int64(v_["T"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("ACX"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ACX"*string(I))
            iv = ix_["VX"*string(I)]
            pbm.A[ig,iv] += Float64(v_["M/H"])
            iv = ix_["VX"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-M/H"])
            iv = ix_["UX"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ACY"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ACY"*string(I))
            iv = ix_["VY"*string(I)]
            pbm.A[ig,iv] += Float64(v_["M/H"])
            iv = ix_["VY"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-M/H"])
            iv = ix_["UY"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("ACZ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"ACZ"*string(I))
            iv = ix_["VZ"*string(I)]
            pbm.A[ig,iv] += Float64(v_["M/H"])
            iv = ix_["VZ"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-M/H"])
            iv = ix_["UZ"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("PSX"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PSX"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            iv = ix_["VX"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("PSY"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PSY"*string(I))
            iv = ix_["Y"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["Y"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            iv = ix_["VY"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("PSZ"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"PSZ"*string(I))
            iv = ix_["Z"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H"])
            iv = ix_["Z"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            iv = ix_["VZ"*string(I)]
            pbm.A[ig,iv] += Float64(-1.0)
            ig,ig_,_ = s2mpj_ii("SC"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"SC"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["T"])
            pbm.gconst[ig_["OX"*string(I)]] = Float64(v_["XT"])
            pbm.gconst[ig_["SC"*string(I)]] = Float64(v_["TOLSQ"])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["Y"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["Y"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["Z"*string(Int64(v_["0"]))]] = 1.0
        pb.xupper[ix_["Z"*string(Int64(v_["0"]))]] = 1.0
        pb.xlower[ix_["VX"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["VX"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["VY"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["VY"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["VZ"*string(Int64(v_["0"]))]] = 0.0
        pb.xupper[ix_["VZ"*string(Int64(v_["0"]))]] = 0.0
        pb.xlower[ix_["VX"*string(Int64(v_["T"]))]] = 0.0
        pb.xupper[ix_["VX"*string(Int64(v_["T"]))]] = 0.0
        pb.xlower[ix_["VY"*string(Int64(v_["T"]))]] = 0.0
        pb.xupper[ix_["VY"*string(Int64(v_["T"]))]] = 0.0
        pb.xlower[ix_["VZ"*string(Int64(v_["T"]))]] = 0.0
        pb.xupper[ix_["VZ"*string(Int64(v_["T"]))]] = 0.0
        for I = Int64(v_["1"]):Int64(v_["T"])
            pb.xlower[ix_["UX"*string(I)]] = v_["-FMAX"]
            pb.xupper[ix_["UX"*string(I)]] = v_["FMAX"]
            pb.xlower[ix_["UY"*string(I)]] = v_["-FMAX"]
            pb.xupper[ix_["UY"*string(I)]] = v_["FMAX"]
            pb.xlower[ix_["UZ"*string(I)]] = v_["-FMAX"]
            pb.xupper[ix_["UZ"*string(I)]] = v_["FMAX"]
            pb.xlower[ix_["X"*string(I)]] = 0.0
            pb.xupper[ix_["X"*string(I)]] = v_["XT"]
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["Y"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["Z"*string(Int64(v_["0"]))]] = Float64(1.0)
        pb.x0[ix_["VX"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["VY"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["VZ"*string(Int64(v_["0"]))]] = Float64(0.0)
        pb.x0[ix_["VX"*string(Int64(v_["T"]))]] = Float64(0.0)
        pb.x0[ix_["VY"*string(Int64(v_["T"]))]] = Float64(0.0)
        pb.x0[ix_["VZ"*string(Int64(v_["T"]))]] = Float64(0.0)
        for I = Int64(v_["1"]):Int64(v_["T"])
            v_["RI"] = Float64(I)
            v_["TI"] = v_["RI"]*v_["H"]
            if haskey(ix_,"X"*string(I))
                pb.x0[ix_["X"*string(I)]] = Float64(v_["TI"])
            else
                pb.y0[findfirst(x->x==ig_["X"*string(I)],pbm.congrps)] = Float64(v_["TI"])
            end
            if haskey(ix_,"VX"*string(I))
                pb.x0[ix_["VX"*string(I)]] = Float64(1.0)
            else
                pb.y0[findfirst(x->x==ig_["VX"*string(I)],pbm.congrps)] = Float64(1.0)
            end
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eERRSIN", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        it,iet_,_ = s2mpj_ii( "eERRCOS", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Z")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["T"])
            ename = "ES"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eERRSIN")
            arrset(ielftype,ie,iet_["eERRSIN"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Y"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            ename = "EC"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eERRCOS")
            arrset(ielftype,ie,iet_["eERRCOS"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "Z"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Z",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gL2",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["T"])
            ig = ig_["OX"*string(I)]
            arrset(pbm.grftype,ig,"gL2")
            ig = ig_["SC"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["ES"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["EC"*string(I)])
            loaset(pbm.grelw,ig,posel, 1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           1.1601050195
# LO SOLTN(50)           26.484181830
# LO SOLTN(100)          44.368110588
# LO SOLTN(500)
# LO SOLTN(1000)
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
        pb.pbclass = "C-SOR2-AN-V-V"
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

    elseif action == "eERRSIN"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SINX = sin(EV_[1])
        COSX = cos(EV_[1])
        ERR = EV_[2]-SINX
        f_   = ERR*ERR
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -2.0*ERR*COSX
            g_[2] = 2.0*ERR
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*(COSX^2+ERR*SINX)
                H_[1,2] = -2.0*COSX
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eERRCOS"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        SINX = sin(EV_[1])
        COSX = cos(EV_[1])
        ERR = EV_[2]-COSX
        f_   = ERR*ERR
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*ERR*SINX
            g_[2] = 2.0*ERR
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*(SINX^2+ERR*COSX)
                H_[1,2] = 2.0*SINX
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gL2"

        GVAR_   = args[1]
        igr_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_= GVAR_*GVAR_
        if nargout>1
            g_ = GVAR_+GVAR_
            if nargout>2
                H_ = zeros(Float64,1,1)
                H_ = 2.0
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

