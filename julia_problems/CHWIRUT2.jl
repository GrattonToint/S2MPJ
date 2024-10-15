function CHWIRUT2(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : CHWIRUT2
#    *********
# 
#    NIST Data fitting problem CHWIRUT2 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = exp[-b1*x]/(b2+b3*x) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Chwirut, D., NIST (197?).  
#      Ultrasonic Reference Block Study. 
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-NOR2-MN-3-54"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "CHWIRUT2"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 54
        v_["N"] = 3
        v_["1"] = 1
        v_["X1"] = 0.5
        v_["X2"] = 1.0
        v_["X3"] = 1.75
        v_["X4"] = 3.75
        v_["X5"] = 5.75
        v_["X6"] = 0.875
        v_["X7"] = 2.25
        v_["X8"] = 3.25
        v_["X9"] = 5.25
        v_["X10"] = 0.75
        v_["X11"] = 1.75
        v_["X12"] = 2.75
        v_["X13"] = 4.75
        v_["X14"] = 0.625
        v_["X15"] = 1.25
        v_["X16"] = 2.25
        v_["X17"] = 4.25
        v_["X18"] = 0.5
        v_["X19"] = 3.0
        v_["X20"] = 0.75
        v_["X21"] = 3.0
        v_["X22"] = 1.5
        v_["X23"] = 6.0
        v_["X24"] = 3.0
        v_["X25"] = 6.0
        v_["X26"] = 1.5
        v_["X27"] = 3.0
        v_["X28"] = 0.5
        v_["X29"] = 2.0
        v_["X30"] = 4.0
        v_["X31"] = 0.75
        v_["X32"] = 2.0
        v_["X33"] = 5.0
        v_["X34"] = 0.75
        v_["X35"] = 2.25
        v_["X36"] = 3.75
        v_["X37"] = 5.75
        v_["X38"] = 3.0
        v_["X39"] = 0.75
        v_["X40"] = 2.5
        v_["X41"] = 4.0
        v_["X42"] = 0.75
        v_["X43"] = 2.5
        v_["X44"] = 4.0
        v_["X45"] = 0.75
        v_["X46"] = 2.5
        v_["X47"] = 4.0
        v_["X48"] = 0.5
        v_["X49"] = 6.0
        v_["X50"] = 3.0
        v_["X51"] = 0.5
        v_["X52"] = 2.75
        v_["X53"] = 0.5
        v_["X54"] = 1.75
        v_["Y1"] = 92.9
        v_["Y2"] = 57.1
        v_["Y3"] = 31.05
        v_["Y4"] = 11.5875
        v_["Y5"] = 8.025
        v_["Y6"] = 63.6
        v_["Y7"] = 21.4
        v_["Y8"] = 14.25
        v_["Y9"] = 8.475
        v_["Y10"] = 63.8
        v_["Y11"] = 26.8
        v_["Y12"] = 16.4625
        v_["Y13"] = 7.125
        v_["Y14"] = 67.3
        v_["Y15"] = 41.0
        v_["Y16"] = 21.15
        v_["Y17"] = 8.175
        v_["Y18"] = 81.50
        v_["Y19"] = 13.12
        v_["Y20"] = 59.9
        v_["Y21"] = 14.62
        v_["Y22"] = 32.9
        v_["Y23"] = 5.44
        v_["Y24"] = 12.56
        v_["Y25"] = 5.44
        v_["Y26"] = 32.0
        v_["Y27"] = 13.95
        v_["Y28"] = 75.8
        v_["Y29"] = 20.0
        v_["Y30"] = 10.42
        v_["Y31"] = 59.5
        v_["Y32"] = 21.67
        v_["Y33"] = 8.55
        v_["Y34"] = 62.0
        v_["Y35"] = 20.2
        v_["Y36"] = 7.76
        v_["Y37"] = 3.75
        v_["Y38"] = 11.81
        v_["Y39"] = 54.7
        v_["Y40"] = 23.7
        v_["Y41"] = 11.55
        v_["Y42"] = 61.3
        v_["Y43"] = 17.7
        v_["Y44"] = 8.74
        v_["Y45"] = 59.2
        v_["Y46"] = 16.3
        v_["Y47"] = 8.62
        v_["Y48"] = 81.0
        v_["Y49"] = 4.87
        v_["Y50"] = 14.62
        v_["Y51"] = 81.7
        v_["Y52"] = 17.17
        v_["Y53"] = 81.3
        v_["Y54"] = 28.9
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
            pb.x0[ix_["B1"]] = Float64(0.1)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(0.1)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(0.01)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(0.01)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(0.02)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(0.02)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE16", iet_)
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
            arrset(pbm.elftype,ie,"eE16")
            arrset(ielftype,ie,iet_["eE16"])
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
        pb.pbclass = "C-NOR2-MN-3-54"
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

    elseif action == "eE16"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        E = exp(-EV_[1]*pbm.elpar[iel_][1])
        EX = E*pbm.elpar[iel_][1]
        EX2 = EX*pbm.elpar[iel_][1]
        V2PV3X = EV_[2]+EV_[3]*pbm.elpar[iel_][1]
        V2PV32 = V2PV3X*V2PV3X
        V2PV33 = V2PV3X*V2PV32
        f_   = E/V2PV3X
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -EX/V2PV3X
            g_[2] = -E/V2PV32
            g_[3] = -EX/V2PV32
            if nargout>2
                H_ = zeros(Float64,3,3)
                H_[1,1] = EX2/V2PV3X
                H_[1,2] = EX/V2PV32
                H_[2,1] = H_[1,2]
                H_[1,3] = EX2/V2PV32
                H_[3,1] = H_[1,3]
                H_[2,2] = 2.0*E/V2PV33
                H_[2,3] = 2.0*EX/V2PV33
                H_[3,2] = H_[2,3]
                H_[3,3] = 2.0*EX2/V2PV33
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

