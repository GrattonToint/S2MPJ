function AIRPORT(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem:
#    ********
# 
#    This problem is concerned with the localisation of airports in Brazil.
#    We consider  m  balls in the real plane, whose centers are the coordinates
#    of some Brazilian  cities and whose  radius were chosen such that the balls are
#    disjoint. The problem is to find one point  (xi, yi) on  each ball, i=1,..,m,
#    such that  SUM(||(xi,yi) - (xj,yj)||)  is  minimum, where the sum involves all
#    the pairs (i,j) such that 1 <= i <= m, 1 <= j <= m and i <> j.
# 
#    For this problem instance, we have m =  42 cities and n = 84 points, 
#    i.e, 42 nonlinear inequalities constraints and 84 variables.
# 
#    Source:
#    Contribution from a LANCELOT user.
# 
#    SIF input : Rodrigo de Barros Nabholz & Maria Aparecida Diniz Ehrhardt
#                November 1994, DMA - IMECC- UNICAMP
#    Adaptation for CUTE: Ph. Toint, November 1994.
# 
#    classification = "C-SQR2-MN-84-42"
# 
#    Problem data
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 6 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "AIRPORT"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 42
        v_["N-1"] = 41
        v_["1"] = 1
        v_["R1"] = 0.09
        v_["R2"] = 0.3
        v_["R3"] = 0.09
        v_["R4"] = 0.45
        v_["R5"] = 0.5
        v_["R6"] = 0.04
        v_["R7"] = 0.1
        v_["R8"] = 0.02
        v_["R9"] = 0.02
        v_["R10"] = 0.07
        v_["R11"] = 0.4
        v_["R12"] = 0.045
        v_["R13"] = 0.05
        v_["R14"] = 0.056
        v_["R15"] = 0.36
        v_["R16"] = 0.08
        v_["R17"] = 0.07
        v_["R18"] = 0.36
        v_["R19"] = 0.67
        v_["R20"] = 0.38
        v_["R21"] = 0.37
        v_["R22"] = 0.05
        v_["R23"] = 0.4
        v_["R24"] = 0.66
        v_["R25"] = 0.05
        v_["R26"] = 0.07
        v_["R27"] = 0.08
        v_["R28"] = 0.3
        v_["R29"] = 0.31
        v_["R30"] = 0.49
        v_["R31"] = 0.09
        v_["R32"] = 0.46
        v_["R33"] = 0.12
        v_["R34"] = 0.07
        v_["R35"] = 0.07
        v_["R36"] = 0.09
        v_["R37"] = 0.05
        v_["R38"] = 0.13
        v_["R39"] = 0.16
        v_["R40"] = 0.46
        v_["R41"] = 0.25
        v_["R42"] = 0.1
        v_["CX1"] = -6.3
        v_["CX2"] = -7.8
        v_["CX3"] = -9.0
        v_["CX4"] = -7.2
        v_["CX5"] = -5.7
        v_["CX6"] = -1.9
        v_["CX7"] = -3.5
        v_["CX8"] = -0.5
        v_["CX9"] = 1.4
        v_["CX10"] = 4.0
        v_["CX11"] = 2.1
        v_["CX12"] = 5.5
        v_["CX13"] = 5.7
        v_["CX14"] = 5.7
        v_["CX15"] = 3.8
        v_["CX16"] = 5.3
        v_["CX17"] = 4.7
        v_["CX18"] = 3.3
        v_["CX19"] = 0.0
        v_["CX20"] = -1.0
        v_["CX21"] = -0.4
        v_["CX22"] = 4.2
        v_["CX23"] = 3.2
        v_["CX24"] = 1.7
        v_["CX25"] = 3.3
        v_["CX26"] = 2.0
        v_["CX27"] = 0.7
        v_["CX28"] = 0.1
        v_["CX29"] = -0.1
        v_["CX30"] = -3.5
        v_["CX31"] = -4.0
        v_["CX32"] = -2.7
        v_["CX33"] = -0.5
        v_["CX34"] = -2.9
        v_["CX35"] = -1.2
        v_["CX36"] = -0.4
        v_["CX37"] = -0.1
        v_["CX38"] = -1.0
        v_["CX39"] = -1.7
        v_["CX40"] = -2.1
        v_["CX41"] = -1.8
        v_["CX42"] = 0.0
        v_["CY1"] = 8.0
        v_["CY2"] = 5.1
        v_["CY3"] = 2.0
        v_["CY4"] = 2.6
        v_["CY5"] = 5.5
        v_["CY6"] = 7.1
        v_["CY7"] = 5.9
        v_["CY8"] = 6.6
        v_["CY9"] = 6.1
        v_["CY10"] = 5.6
        v_["CY11"] = 4.9
        v_["CY12"] = 4.7
        v_["CY13"] = 4.3
        v_["CY14"] = 3.6
        v_["CY15"] = 4.1
        v_["CY16"] = 3.0
        v_["CY17"] = 2.4
        v_["CY18"] = 3.0
        v_["CY19"] = 4.7
        v_["CY20"] = 3.4
        v_["CY21"] = 2.3
        v_["CY22"] = 1.5
        v_["CY23"] = 0.5
        v_["CY24"] = -1.7
        v_["CY25"] = -2.0
        v_["CY26"] = -3.1
        v_["CY27"] = -3.5
        v_["CY28"] = -2.4
        v_["CY29"] = -1.3
        v_["CY30"] = 0.0
        v_["CY31"] = -1.7
        v_["CY32"] = -2.1
        v_["CY33"] = -0.4
        v_["CY34"] = -2.9
        v_["CY35"] = -3.4
        v_["CY36"] = -4.3
        v_["CY37"] = -5.2
        v_["CY38"] = -6.5
        v_["CY39"] = -7.5
        v_["CY40"] = -6.4
        v_["CY41"] = -5.1
        v_["CY42"] = 0.0
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
            iv,ix_,_ = s2mpj_ii("Y"*string(I),ix_)
            arrset(pb.xnames,iv,"Y"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = I+v_["1"]
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig,ig_,_ = s2mpj_ii("OBJ1"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["X"*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["X"*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
                ig,ig_,_ = s2mpj_ii("OBJ2"*string(I)*","*string(J),ig_)
                arrset(gtype,ig,"<>")
                iv = ix_["Y"*string(I)]
                pbm.A[ig,iv] += Float64(1.0)
                iv = ix_["Y"*string(J)]
                pbm.A[ig,iv] += Float64(-1.0)
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig,ig_,_ = s2mpj_ii("CONS"*string(I),ig_)
            arrset(gtype,ig,"<=")
            arrset(pb.cnames,ig,"CONS"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["N"])
            pbm.gconst[ig_["CONS"*string(I)]] = Float64(v_["R"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        for I = Int64(v_["1"]):Int64(v_["N"])
            pb.xlower[ix_["X"*string(I)]] = -10
            pb.xupper[ix_["X"*string(I)]] = 10
            pb.xlower[ix_["Y"*string(I)]] = -10
            pb.xupper[ix_["Y"*string(I)]] = 10
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eDIFSQR", iet_)
        loaset(elftv,it,1,"V")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"W")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "A"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDIFSQR")
            arrset(ielftype,ie,iet_["eDIFSQR"])
            vname = "X"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="W",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CX"*string(I)]))
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ename = "B"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eDIFSQR")
            arrset(ielftype,ie,iet_["eDIFSQR"])
            vname = "Y"*string(I)
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="W",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["CY"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = Dict{String,Int}()
        it,igt_,_ = s2mpj_ii("gSQUARE",igt_)
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N-1"])
            v_["I+1"] = I+v_["1"]
            for J = Int64(v_["I+1"]):Int64(v_["N"])
                ig = ig_["OBJ1"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gSQUARE")
                ig = ig_["OBJ2"*string(I)*","*string(J)]
                arrset(pbm.grftype,ig,"gSQUARE")
            end
        end
        for I = Int64(v_["1"]):Int64(v_["N"])
            ig = ig_["CONS"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["A"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,Float64(1.0))
            posel = posel+1
            loaset(pbm.grelt,ig,posel,ie_["B"*string(I)])
            loaset(pbm.grelw,ig,posel,Float64(1.0))
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = .0
#    Solution
# LO SOLTN              47952.695811
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.cupper[1:pb.nle] = zeros(Float64,pb.nle)
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-SQR2-MN-84-42"
        pb.x0          = zeros(Float64,pb.n)
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

    elseif action == "eDIFSQR"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DIF = EV_[1]-pbm.elpar[iel_][1]
        f_   = DIF*DIF
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 2.0*DIF
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

    #%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    elseif action == "gSQUARE"

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

