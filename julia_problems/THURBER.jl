function THURBER(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : THURBER
#    *********
# 
#    NIST Data fitting problem THURBER given as an inconsistent set of
#    nonlinear equations.
# 
#    Fit: y = (b1 + b2*x + b3*x**2 + b4*x**3) / 
#             (1 + b5*x + b6*x**2 + b7*x**3) + e
# 
#    Source:  Problem from the NIST nonlinear regression test set
#      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
# 
#    Reference: Thurber, R., NIST (197?).  
#      Semiconductor electron mobility modeling.
# 
#    SIF input: Nick Gould and Tyrone Rees, Oct 2015
# 
#    classification = "C-NOR2-MN-7-37"
# 
#    Number of data values
# 
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Translated to Julia by S2MPJ version 7 X 2024
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "THURBER"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["M"] = 37
        v_["N"] = 7
        v_["1"] = 1
        v_["X1"] = -3.067
        v_["X2"] = -2.981
        v_["X3"] = -2.921
        v_["X4"] = -2.912
        v_["X5"] = -2.840
        v_["X6"] = -2.797
        v_["X7"] = -2.702
        v_["X8"] = -2.699
        v_["X9"] = -2.633
        v_["X10"] = -2.481
        v_["X11"] = -2.363
        v_["X12"] = -2.322
        v_["X13"] = -1.501
        v_["X14"] = -1.460
        v_["X15"] = -1.274
        v_["X16"] = -1.212
        v_["X17"] = -1.100
        v_["X18"] = -1.046
        v_["X19"] = -0.915
        v_["X20"] = -0.714
        v_["X21"] = -0.566
        v_["X22"] = -0.545
        v_["X23"] = -0.400
        v_["X24"] = -0.309
        v_["X25"] = -0.109
        v_["X26"] = -0.103
        v_["X27"] = 0.010
        v_["X28"] = 0.119
        v_["X29"] = 0.377
        v_["X30"] = 0.790
        v_["X31"] = 0.963
        v_["X32"] = 1.006
        v_["X33"] = 1.115
        v_["X34"] = 1.572
        v_["X35"] = 1.841
        v_["X36"] = 2.047
        v_["X37"] = 2.200
        v_["Y1"] = 80.574
        v_["Y2"] = 84.248
        v_["Y3"] = 87.264
        v_["Y4"] = 87.195
        v_["Y5"] = 89.076
        v_["Y6"] = 89.608
        v_["Y7"] = 89.868
        v_["Y8"] = 90.101
        v_["Y9"] = 92.405
        v_["Y10"] = 95.854
        v_["Y11"] = 100.696
        v_["Y12"] = 101.060
        v_["Y13"] = 401.672
        v_["Y14"] = 390.724
        v_["Y15"] = 567.534
        v_["Y16"] = 635.316
        v_["Y17"] = 733.054
        v_["Y18"] = 759.087
        v_["Y19"] = 894.206
        v_["Y20"] = 990.785
        v_["Y21"] = 1090.109
        v_["Y22"] = 1080.914
        v_["Y23"] = 1122.643
        v_["Y24"] = 1178.351
        v_["Y25"] = 1260.531
        v_["Y26"] = 1273.514
        v_["Y27"] = 1288.339
        v_["Y28"] = 1327.543
        v_["Y29"] = 1353.863
        v_["Y30"] = 1414.509
        v_["Y31"] = 1425.208
        v_["Y32"] = 1421.384
        v_["Y33"] = 1442.962
        v_["Y34"] = 1464.350
        v_["Y35"] = 1468.705
        v_["Y36"] = 1447.894
        v_["Y37"] = 1457.628
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
            pb.x0[ix_["B1"]] = Float64(1000.0)
        else
            pb.y0[findfirst(x->x==ig_["B1"],pbm.congrps)] = Float64(1000.0)
        end
        if haskey(ix_,"B2")
            pb.x0[ix_["B2"]] = Float64(1000.0)
        else
            pb.y0[findfirst(x->x==ig_["B2"],pbm.congrps)] = Float64(1000.0)
        end
        if haskey(ix_,"B3")
            pb.x0[ix_["B3"]] = Float64(400.0)
        else
            pb.y0[findfirst(x->x==ig_["B3"],pbm.congrps)] = Float64(400.0)
        end
        if haskey(ix_,"B4")
            pb.x0[ix_["B4"]] = Float64(40.0)
        else
            pb.y0[findfirst(x->x==ig_["B4"],pbm.congrps)] = Float64(40.0)
        end
        if haskey(ix_,"B5")
            pb.x0[ix_["B5"]] = Float64(0.7)
        else
            pb.y0[findfirst(x->x==ig_["B5"],pbm.congrps)] = Float64(0.7)
        end
        if haskey(ix_,"B6")
            pb.x0[ix_["B6"]] = Float64(0.3)
        else
            pb.y0[findfirst(x->x==ig_["B6"],pbm.congrps)] = Float64(0.3)
        end
        if haskey(ix_,"B7")
            pb.x0[ix_["B7"]] = Float64(0.03)
        else
            pb.y0[findfirst(x->x==ig_["B7"],pbm.congrps)] = Float64(0.03)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eE19", iet_)
        loaset(elftv,it,1,"V1")
        loaset(elftv,it,2,"V2")
        loaset(elftv,it,3,"V3")
        loaset(elftv,it,4,"V4")
        loaset(elftv,it,5,"V5")
        loaset(elftv,it,6,"V6")
        loaset(elftv,it,7,"V7")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"X")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["M"])
            ename = "E"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eE19")
            arrset(ielftype,ie,iet_["eE19"])
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
            vname = "B4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B6"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V6",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "B7"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="V7",elftv[ielftype[ie]])
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
        pb.pbclass = "C-NOR2-MN-7-37"
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

    elseif action == "eE19"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        X2 = pbm.elpar[iel_][1]*pbm.elpar[iel_][1]
        X3 = X2*pbm.elpar[iel_][1]
        X4 = X3*pbm.elpar[iel_][1]
        X5 = X4*pbm.elpar[iel_][1]
        X6 = X5*pbm.elpar[iel_][1]
        T = EV_[1]+EV_[2]*pbm.elpar[iel_][1]+EV_[3]*X2+EV_[4]*X3
        D = 1.0e0+EV_[5]*pbm.elpar[iel_][1]+EV_[6]*X2+EV_[7]*X3
        D2 = D*D
        TD3 = 0.5e0*D2*D
        f_   = T/D
        if nargout>1
            dim = try length(IV_) catch; length(EV_) end
            g_  = zeros(Float64,dim)
            g_[1] = 1.0e0/D
            g_[2] = pbm.elpar[iel_][1]/D
            g_[3] = X2/D
            g_[4] = X3/D
            g_[5] = -pbm.elpar[iel_][1]*T/D2
            g_[6] = -X2*T/D2
            g_[7] = -X3*T/D2
            if nargout>2
                H_ = zeros(Float64,7,7)
                H_[1,5] = -pbm.elpar[iel_][1]/D2
                H_[5,1] = H_[1,5]
                H_[1,6] = -X2/D2
                H_[6,1] = H_[1,6]
                H_[1,7] = -X3/D2
                H_[7,1] = H_[1,7]
                H_[2,5] = -X2/D2
                H_[5,2] = H_[2,5]
                H_[2,6] = -X3/D2
                H_[6,2] = H_[2,6]
                H_[2,7] = -X4/D2
                H_[7,2] = H_[2,7]
                H_[3,5] = -X3/D2
                H_[5,3] = H_[3,5]
                H_[3,6] = -X4/D2
                H_[6,3] = H_[3,6]
                H_[3,7] = -X5/D2
                H_[7,3] = H_[3,7]
                H_[4,5] = -X4/D2
                H_[5,4] = H_[4,5]
                H_[4,6] = -X5/D2
                H_[6,4] = H_[4,6]
                H_[4,7] = -X6/D2
                H_[7,4] = H_[4,7]
                H_[5,5] = X2*T/TD3
                H_[5,6] = X3*T/TD3
                H_[6,5] = H_[5,6]
                H_[5,7] = X4*T/TD3
                H_[7,5] = H_[5,7]
                H_[6,6] = X4*T/TD3
                H_[6,7] = X5*T/TD3
                H_[7,6] = H_[6,7]
                H_[7,7] = X6*T/TD3
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

