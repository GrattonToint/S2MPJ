function ROSZMAN1(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : ROSZMAN1
#    *********
# 
#    NIST Data fitting problem ROSZMAN1 given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y =  b1 - b2*x - arctan[b3/(x-b4)]/pi + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#   Reference: Roszman, L., NIST (1979).  
#     Quantum Defects for Sulfur I Atom.
# 
#    classification = "C-NOR2-MN-4-25"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "ROSZMAN1"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 25
        v_["N"] = 4
        v_["1"] = 1
        v_["X1"] = -4868.68
        v_["X2"] = -4868.09
        v_["X3"] = -4867.41
        v_["X4"] = -3375.19
        v_["X5"] = -3373.14
        v_["X6"] = -3372.03
        v_["X7"] = -2473.74
        v_["X8"] = -2472.35
        v_["X9"] = -2469.45
        v_["X10"] = -1894.65
        v_["X11"] = -1893.40
        v_["X12"] = -1497.24
        v_["X13"] = -1495.85
        v_["X14"] = -1493.41
        v_["X15"] = -1208.68
        v_["X16"] = -1206.18
        v_["X17"] = -1206.04
        v_["X18"] = -997.92
        v_["X19"] = -996.61
        v_["X20"] = -996.31
        v_["X21"] = -834.94
        v_["X22"] = -834.66
        v_["X23"] = -710.03
        v_["X24"] = -530.16
        v_["X25"] = -464.17
        v_["Y1"] = 0.252429
        v_["Y2"] = 0.252141
        v_["Y3"] = 0.251809
        v_["Y4"] = 0.297989
        v_["Y5"] = 0.296257
        v_["Y6"] = 0.295319
        v_["Y7"] = 0.339603
        v_["Y8"] = 0.337731
        v_["Y9"] = 0.333820
        v_["Y10"] = 0.389510
        v_["Y11"] = 0.386998
        v_["Y12"] = 0.438864
        v_["Y13"] = 0.434887
        v_["Y14"] = 0.427893
        v_["Y15"] = 0.471568
        v_["Y16"] = 0.461699
        v_["Y17"] = 0.461144
        v_["Y18"] = 0.513532
        v_["Y19"] = 0.506641
        v_["Y20"] = 0.505062
        v_["Y21"] = 0.535648
        v_["Y22"] = 0.533726
        v_["Y23"] = 0.568064
        v_["Y24"] = 0.612886
        v_["Y25"] = 0.624169
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
            iv = ix_["B1"]
            pbm.A[ig,iv] += Float64(1.0)
            v_["-X"] = -1.0*v_["X"*string(I)]
            iv = ix_["B2"]
            pbm.A[ig,iv] += Float64(v_["-X"])
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
            pb.x0[ix_["B2"]] = Float64(-0.00001)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(-0.00001)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(1000.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(1000.0)
        end
        if haskey(ix_,"B4")
            pb.x0[ix_["B4"]] = Float64(-100.0)
        else
            pb.y0[findfirst(x->x==ig_["B4"],pbm.congrps)] = Float64(-100.0)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE7", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE7")
            arrset(ielftype,ie,iet_["eE7"])
            vname = "B3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V2",elftv[ielftype[ie]])
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
        Asave = pbm.A[1:ngrp, 1:pb.n]
        pbm.A = Asave
        pbm.H = spzeros(Float64,0,0)
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-NOR2-MN-4-25"
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

    elseif action == "e_globs"

        pbm = args[1]
        arrset(pbm.efpar,1,4.0*atan(1.0e0))
        return pbm

    elseif action == "eE7"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        V12 = EV_[1]*EV_[1]
        V13 = EV_[1]*V12
        V2MX = EV_[2]-pbm.elpar[iel_][1]
        V2MX2 = V2MX*V2MX
        V2MX3 = V2MX*V2MX2
        R = V12/V2MX2+1.0
        PIR = pbm.efpar[1]*R
        PIR2 = PIR*R
        f_   = -atan(EV_[1]/V2MX)/pbm.efpar[1]
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = -1.0/(PIR*V2MX)
            g_[2] = EV_[1]/(PIR*V2MX2)
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,1] = 2.0*EV_[1]/(PIR2*V2MX3)
                H_[1,2] = 1.0/(PIR*V2MX2)-2.0*V12/(PIR2*V2MX^4)
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*V13/(PIR2*V2MX^5)-2.0*EV_[1]/(PIR*V2MX3)
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
            pbm.has_globs = [1,0]
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

