function PALMER1BNE(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : PALMER1BNE
#    *********
# 
#    A nonlinear least squares problem with bounds
#    arising from chemical kinetics.
# 
#    model: H-N=N=N TZVP+MP2
#    fitting Y to A2 X**2 + A4 X**4
#                 + B / ( C + X**2 ), B, C nonnegative.
# 
#    Source:
#    M. Palmer, Edinburgh, private communication.
# 
#    SIF input: Nick Gould, 1990.
#    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
# 
#    classification = "C-NOR2-RN-4-35"
# 
#    Number of data points
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "PALMER1BNE"

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
        v_["1"] = 1
        v_["X1"] = -1.788963
        v_["X2"] = -1.745329
        v_["X3"] = -1.658063
        v_["X4"] = -1.570796
        v_["X5"] = -1.483530
        v_["X6"] = -1.396263
        v_["X7"] = -1.308997
        v_["X8"] = -1.218612
        v_["X9"] = -1.134464
        v_["X10"] = -1.047198
        v_["X11"] = -0.872665
        v_["X12"] = -0.698132
        v_["X13"] = -0.523599
        v_["X14"] = -0.349066
        v_["X15"] = -0.174533
        v_["X16"] = 0.0000000
        v_["X17"] = 1.788963
        v_["X18"] = 1.745329
        v_["X19"] = 1.658063
        v_["X20"] = 1.570796
        v_["X21"] = 1.483530
        v_["X22"] = 1.396263
        v_["X23"] = 1.308997
        v_["X24"] = 1.218612
        v_["X25"] = 1.134464
        v_["X26"] = 1.047198
        v_["X27"] = 0.872665
        v_["X28"] = 0.698132
        v_["X29"] = 0.523599
        v_["X30"] = 0.349066
        v_["X31"] = 0.174533
        v_["X32"] = -1.8762289
        v_["X33"] = -1.8325957
        v_["X34"] = 1.8762289
        v_["X35"] = 1.8325957
        v_["Y1"] = 78.596218
        v_["Y2"] = 65.77963
        v_["Y3"] = 43.96947
        v_["Y4"] = 27.038816
        v_["Y5"] = 14.6126
        v_["Y6"] = 6.2614
        v_["Y7"] = 1.538330
        v_["Y8"] = 0.000000
        v_["Y9"] = 1.188045
        v_["Y10"] = 4.6841
        v_["Y11"] = 16.9321
        v_["Y12"] = 33.6988
        v_["Y13"] = 52.3664
        v_["Y14"] = 70.1630
        v_["Y15"] = 83.4221
        v_["Y16"] = 88.3995
        v_["Y17"] = 78.596218
        v_["Y18"] = 65.77963
        v_["Y19"] = 43.96947
        v_["Y20"] = 27.038816
        v_["Y21"] = 14.6126
        v_["Y22"] = 6.2614
        v_["Y23"] = 1.538330
        v_["Y24"] = 0.000000
        v_["Y25"] = 1.188045
        v_["Y26"] = 4.6841
        v_["Y27"] = 16.9321
        v_["Y28"] = 33.6988
        v_["Y29"] = 52.3664
        v_["Y30"] = 70.1630
        v_["Y31"] = 83.4221
        v_["Y32"] = 108.18086
        v_["Y33"] = 92.733676
        v_["Y34"] = 108.18086
        v_["Y35"] = 92.733676
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        iv,ix_,_ = s2mpj_ii("A2",ix_)
        arrset(pb.xnames,iv,"A2")
        iv,ix_,_ = s2mpj_ii("A4",ix_)
        arrset(pb.xnames,iv,"A4")
        iv,ix_,_ = s2mpj_ii("B",ix_)
        arrset(pb.xnames,iv,"B")
        iv,ix_,_ = s2mpj_ii("C",ix_)
        arrset(pb.xnames,iv,"C")
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype    = String[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            v_["XQUART"] = v_["XSQR"]*v_["XSQR"]
            ig,ig_,_ = s2mpj_ii("O"*string(I),ig_)
            arrset(gtype,ig,"==")
            arrset(pb.cnames,ig,"O"*string(I))
            iv = ix_["A2"]
            pbm.A[ig,iv] += Float64(v_["XSQR"])
            iv = ix_["A4"]
            pbm.A[ig,iv] += Float64(v_["XQUART"])
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
            pbm.gconst[ig_["O"*string(I)]] = Float64(v_["Y"*string(I)])
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["A2"]] = -Inf
        pb.xupper[ix_["A2"]] = +Inf
        pb.xlower[ix_["A4"]] = -Inf
        pb.xupper[ix_["A4"]] = +Inf
        pb.xlower[ix_["B"]] = 0.00001
        pb.xlower[ix_["C"]] = 0.00001
        #%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = fill(Float64(1.0),pb.n)
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eQUOT", iet_)
        loaset(elftv,it,1,"B")
        loaset(elftv,it,2,"C")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"XSQR")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            v_["XSQR"] = v_["X"*string(I)]*v_["X"*string(I)]
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eQUOT")
            arrset(ielftype,ie,iet_["eQUOT"])
            vname = "B"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="B",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "C"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,Float64(1.0))
            posev = findfirst(x->x=="C",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="XSQR",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["XSQR"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        for I = Int64(v_["1"]):Int64(v_["M"])
            ig = ig_["O"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["E"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Least square problems are bounded below by zero
# LO PALMER1B               0.0
#    Solution
# LO SOLTN               3.44734948
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
        pb.pbclass = "C-NOR2-RN-4-35"
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

    elseif action == "eQUOT"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        DENOM = 1.0/(EV_[2]+pbm.elpar[iel_][1])
        f_   = EV_[1]*DENOM
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = DENOM
            g_[2] = -EV_[1]*DENOM*DENOM
            if nargout>2
                H_ = zeros(Float64,2,2)
                H_[1,2] = -DENOM*DENOM
                H_[2,1] = H_[1,2]
                H_[2,2] = 2.0*EV_[1]*DENOM^3
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

