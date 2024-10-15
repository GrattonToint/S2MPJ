function HAGER4(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HAGER4
#    *********
# 
#    A nonlinear optimal control problem, by W. Hager.
# 
#    NOTE: The solution for x given in the article below by Hager has
#    a typo. On the interval [1/2, 1], x(t) = (exp(2t) + exp(t))/d. In
#    other words, the minus sign in the article should be a plus sign.
# 
#    Source: problem P4 in
#    W.W. Hager,
#    "Multiplier Methods for Nonlinear Optimal Control",
#    SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.
# 
#    SIF input: Ph. Toint, April 1991.
# 
#    classification = "C-OLR2-AN-V-V"
# 
#    Number of discretized points in [0,1]
# 
#       Alternative values for the SIF file parameters:
# IE N                   10             $-PARAMETER     original value
# IE N                   50             $-PARAMETER
# IE N                   100            $-PARAMETER
# IE N                   500            $-PARAMETER
# IE N                   1000           $-PARAMETER
# IE N                   2500           $-PARAMETER
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HAGER4"

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
# IE N                   5000           $-PARAMETER
        v_["1/H"] = Float64(v_["N"])
        v_["H"] = 1.0/v_["1/H"]
        v_["H/2"] = 0.5*v_["H"]
        v_["1/H-1"] = -1.0+v_["1/H"]
        v_["-1/H"] = -1.0*v_["1/H"]
        v_["1/HSQ"] = v_["1/H"]*v_["1/H"]
        v_["1/2HSQ"] = 0.5*v_["1/HSQ"]
        v_["0"] = 0
        v_["1"] = 1
        for I = Int64(v_["0"]):Int64(v_["N"])
            v_["RI"] = Float64(I)
            v_["T"*string(I)] = v_["RI"]*v_["H"]
            v_["-2TI"] = -2.0*v_["T"*string(I)]
            v_["Z"*string(I)] = exp(v_["-2TI"])
        end
        for I = Int64(v_["0"]):Int64(v_["1"])
            v_["A"*string(I)] = -0.5*v_["Z"*string(I)]
            v_["TI+1/2"] = 0.5+v_["T"*string(I)]
            v_["B"*string(I)] = v_["A"*string(I)]*v_["TI+1/2"]
            v_["TISQ"] = v_["T"*string(I)]*v_["T"*string(I)]
            v_["TIETC"] = v_["TISQ"]+v_["TI+1/2"]
            v_["C"*string(I)] = v_["A"*string(I)]*v_["TIETC"]
        end
        v_["DA"] = v_["A"*string(Int64(v_["1"]))]-v_["A"*string(Int64(v_["0"]))]
        v_["SCDA"] = 0.5*v_["DA"]
        v_["DB"] = v_["B"*string(Int64(v_["1"]))]-v_["B"*string(Int64(v_["0"]))]
        v_["SCDB"] = v_["DB"]*v_["1/H"]
        v_["DC"] = v_["C"*string(Int64(v_["1"]))]-v_["C"*string(Int64(v_["0"]))]
        v_["SCDC"] = v_["DC"]*v_["1/2HSQ"]
        v_["E"] = exp(1.0)
        v_["3E"] = 3.0*v_["E"]
        v_["1+3E"] = 1.0+v_["3E"]
        v_["1-E"] = 1.0-v_["E"]
        v_["2-2E"] = 2.0*v_["1-E"]
        v_["XX0"] = v_["1+3E"]/v_["2-2E"]
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["0"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("U"*string(I),ix_)
            arrset(pb.xnames,iv,"U"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ig,ig_,_ = s2mpj_ii("S"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"S"*string(I))
            iv = ix_["X"*string(I)]
            pbm.A[ig,iv] += Float64(v_["1/H-1"])
            iv = ix_["X"*string(Int64(v_["I-1"]))]
            pbm.A[ig,iv] += Float64(v_["-1/H"])
            v_["ETI"] = exp(v_["T"*string(I)])
            v_["-ETI"] = -1.0*v_["ETI"]
            iv = ix_["U"*string(I)]
            pbm.A[ig,iv] += Float64(v_["-ETI"])
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
        pb.xlower[ix_["X"*string(Int64(v_["0"]))]] = v_["XX0"]
        pb.xupper[ix_["X"*string(Int64(v_["0"]))]] = v_["XX0"]
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xupper[ix_["U"*string(I)]] = 1.0
        end
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        pb.x0[ix_["X"*string(Int64(v_["0"]))]] = Float64(v_["XX0"])
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eELT", iet_)
        loaset(elftv,it,1,"X")
        loaset(elftv,it,2,"Y")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"D")
        loaset(elftp,it,2,"E")
        loaset(elftp,it,3,"F")
        it,iet_,_ = s2mpj_ii( "eSQ", iet_)
        loaset(elftv,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            v_["I-1"] = -1+I
            ename = "EL"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eELT")
            arrset(ielftype,ie,iet_["eELT"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X"*string(Int64(v_["I-1"]))
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="Y",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            v_["DD"] = v_["SCDA"]*v_["Z"*string(Int64(v_["I-1"]))]
            v_["EE"] = v_["SCDB"]*v_["Z"*string(Int64(v_["I-1"]))]
            v_["FF"] = v_["SCDC"]*v_["Z"*string(Int64(v_["I-1"]))]
            posep = findfirst(x->x=="D",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["DD"]))
            posep = findfirst(x->x=="E",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["EE"]))
            posep = findfirst(x->x=="F",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["FF"]))
            ename = "U"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eSQ")
            arrset(ielftype,ie,iet_["eSQ"])
            vname = "U"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["OBJ"]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["EL"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["U"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(v_["H/2"]))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN(10)           2.833914199
# LO SOLTN(50)           2.799810928
# LO SOLTN(100)          2.796761851
# LO SOLTN(500)          2.794513229
# LO SOLTN(1000)         2.794244187
# LO SOLTN(5000)         ???
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
        pb.pbclass = "C-OLR2-AN-V-V"
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

    elseif action == "eELT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        f_    = (
              pbm.elpar[iel_][1]*EV_[1]*EV_[1]+pbm.elpar[iel_][2]*EV_[1]*(EV_[2]-EV_[1])+pbm.elpar[iel_][3]*(EV_[2]-EV_[1])^2)
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = (2.0*pbm.elpar[iel_][1]*EV_[1]+pbm.elpar[iel_][2]*(EV_[2]-2.0*EV_[1])-
                 2.0*pbm.elpar[iel_][3]*(EV_[2]-EV_[1]))
            g_[2] = pbm.elpar[iel_][2]*EV_[1]+2.0*pbm.elpar[iel_][3]*(EV_[2]-EV_[1])
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*(pbm.elpar[iel_][1]-pbm.elpar[iel_][2]+pbm.elpar[iel_][3])
                H_[1,2] = pbm.elpar[iel_][2]-2.0*pbm.elpar[iel_][3]
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*pbm.elpar[iel_][3]
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

