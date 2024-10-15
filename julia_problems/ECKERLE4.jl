function ECKERLE4(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ECKERLE4
#    *********
# 
#    NIST Data fitting problem ECKERLE4 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = (b1/b2) * exp[-0.5*((x-b3)/b2)**2] + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Eckerle, K., NIST (197?).  
#      Circular Interference Transmittance Study.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-NOR2-MN-3-35"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ECKERLE4"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 35
        v_["N"] = 3
        v_["1"] = 1
        v_["X1"] = 400.000000
        v_["X2"] = 405.000000
        v_["X3"] = 410.000000
        v_["X4"] = 415.000000
        v_["X5"] = 420.000000
        v_["X6"] = 425.000000
        v_["X7"] = 430.000000
        v_["X8"] = 435.000000
        v_["X9"] = 436.500000
        v_["X10"] = 438.000000
        v_["X11"] = 439.500000
        v_["X12"] = 441.000000
        v_["X13"] = 442.500000
        v_["X14"] = 444.000000
        v_["X15"] = 445.500000
        v_["X16"] = 447.000000
        v_["X17"] = 448.500000
        v_["X18"] = 450.000000
        v_["X19"] = 451.500000
        v_["X20"] = 453.000000
        v_["X21"] = 454.500000
        v_["X22"] = 456.000000
        v_["X23"] = 457.500000
        v_["X24"] = 459.000000
        v_["X25"] = 460.500000
        v_["X26"] = 462.000000
        v_["X27"] = 463.500000
        v_["X28"] = 465.000000
        v_["X29"] = 470.000000
        v_["X30"] = 475.000000
        v_["X31"] = 480.000000
        v_["X32"] = 485.000000
        v_["X33"] = 490.000000
        v_["X34"] = 495.000000
        v_["X35"] = 500.000000
        v_["Y1"] = 0.0001575
        v_["Y2"] = 0.0001699
        v_["Y3"] = 0.0002350
        v_["Y4"] = 0.0003102
        v_["Y5"] = 0.0004917
        v_["Y6"] = 0.0008710
        v_["Y7"] = 0.0017418
        v_["Y8"] = 0.0046400
        v_["Y9"] = 0.0065895
        v_["Y10"] = 0.0097302
        v_["Y11"] = 0.0149002
        v_["Y12"] = 0.0237310
        v_["Y13"] = 0.0401683
        v_["Y14"] = 0.0712559
        v_["Y15"] = 0.1264458
        v_["Y16"] = 0.2073413
        v_["Y17"] = 0.2902366
        v_["Y18"] = 0.3445623
        v_["Y19"] = 0.3698049
        v_["Y20"] = 0.3668534
        v_["Y21"] = 0.3106727
        v_["Y22"] = 0.2078154
        v_["Y23"] = 0.1164354
        v_["Y24"] = 0.0616764
        v_["Y25"] = 0.0337200
        v_["Y26"] = 0.0194023
        v_["Y27"] = 0.0117831
        v_["Y28"] = 0.0074357
        v_["Y29"] = 0.0022732
        v_["Y30"] = 0.0008800
        v_["Y31"] = 0.0004579
        v_["Y32"] = 0.0002345
        v_["Y33"] = 0.0001586
        v_["Y34"] = 0.0001143
        v_["Y35"] = 0.0000710
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("B"*string(I),ix_)
            arrset(pb.xnames,iv,"B"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig,ig_,_ = s2mpj_ii("F"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"F"*string(I))
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
        for I = Int64(v_["1"]):Int64(v_["M"])
            pbm.gconst[ig_["F"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1*fill(Inf,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower = zeros(Float64,pb.n)
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"B1")
            pb.x0[ix_["B1"]] = Float64(1.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(1.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(10.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(10.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(500.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(500.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE")
            arrset(ielftype,ie,iet_["eE"])
            vname = "B1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="X",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["X"*string(I)]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["F"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
        pb.objlower = 0.0
#    Solution
# LO SOLTN               
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pb.cupper[pb.nle+1:pb.nle+pb.neq] = zeros(Float64,pb.neq)
        pbm.A = spzeros(Float64,0,0)
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-MN-3-35"
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

    elseif action == "eE"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V3MX = EV_[3]-pbm.elpar[iel_][1]
        TV3MX = 2.0*V3MX
        V3MX2 = V3MX^2
        V22 = EV_[2]^2
        V23 = EV_[2]*V22
        V24 = V22*V22
        V25 = V22*V23
        V26 = V23*V23
        V27 = V23*V24
        E = exp(-0.5*V3MX2/V22)
        V1E = EV_[1]*E
        DIFF = V3MX2/V24-1.0/V22
        f_   = EV_[1]*E/EV_[2]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = E/EV_[2]
            g_[2] = V1E*DIFF
            g_[3] = -V1E*V3MX/V23
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,2] = E*DIFF
                H_[2,1] = H_[1,2]
                H_[1,3] = -0.5*E*TV3MX/V23
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*V1E/V23-5.0*V1E*V3MX2/V25+V1E*V3MX^4/V27
                H_[2,3] = 1.5*V1E*TV3MX/V24-0.5*V1E*TV3MX*V3MX2/V26
                H_[3,2] = H_[2,3]
                H_[3,3] = 0.5*V1E*V3MX2/V25-V1E/V23
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

