function HS85(action::String,args::Union{PBM,Int,Float64,Vector{Int},Vector{Float64}}...)
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# 
#    Problem : HS85
#    *********
# 
#    The problem is to optimize the net profit of an hypothetical
#    wood-pulp plant. The constraints include the usual material
#    and energy balances as well as several empirical equations.
# 
#    Source: problem 85 in
#    W. Hock and K. Schittkowski,
#    "Test examples for nonlinear programming codes",
#    Lectures Notes in Economics and Mathematical Systems 187, Springer
#    Verlag, Heidelberg, 1981.
# 
#    SIF input: Nick Gould, September 1991.
#      SAVEs removed December 3rd 2014
#      Julia coding : Cunxin Huang, 2025.
# 
#    classification = "C-COOI2-MN-5-21"
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = "HS85"

    if action == "setup"
        pb           = PB(name)
        pbm          = PBM(name)
        nargin       = length(args)
        pbm.call     = getfield( Main, Symbol( name ) )

        #%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = Dict{String,Float64}();
        ix_ = Dict{String,Int}();
        ig_ = Dict{String,Int}();
        v_["N"] = 5
        v_["A2"] = 17.505
        v_["A3"] = 11.275
        v_["A4"] = 214.228
        v_["A5"] = 7.458
        v_["A6"] = 0.961
        v_["A7"] = 1.612
        v_["A8"] = 0.146
        v_["A9"] = 107.99
        v_["A10"] = 922.693
        v_["A11"] = 926.832
        v_["A12"] = 18.766
        v_["A13"] = 1072.163
        v_["A14"] = 8961.448
        v_["A15"] = 0.063
        v_["A16"] = 71084.33
        v_["A17"] = 2802713.0
        v_["B2"] = 1053.6667
        v_["B3"] = 35.03
        v_["B4"] = 665.585
        v_["B5"] = 584.463
        v_["B6"] = 265.916
        v_["B7"] = 7.046
        v_["B8"] = 0.222
        v_["B9"] = 273.366
        v_["B10"] = 1286.105
        v_["B11"] = 1444.046
        v_["B12"] = 537.141
        v_["B13"] = 3247.039
        v_["B14"] = 26844.086
        v_["B15"] = 0.386
        v_["B16"] = 140000.0
        v_["B17"] = 12146108.0
        v_["1"] = 1
        v_["2"] = 2
        v_["17"] = 17
        v_["19"] = 19
        v_["20"] = 20
        #%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xscale = Float64[]
        intvars = Int64[]
        binvars = Int64[]
        irA   = Int64[]
        icA   = Int64[]
        valA  = Float64[]
        for I = Int64(v_["1"]):Int64(v_["N"])
            iv,ix_,_ = s2mpj_ii("X"*string(I),ix_)
            arrset(pb.xnames,iv,"X"*string(I))
        end
        #%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        gtype = String[]
        ig,ig_,_ = s2mpj_ii("OBJ",ig_)
        arrset(gtype,ig,"<>")
        ig,ig_,_ = s2mpj_ii("CON0",ig_)
        arrset(gtype,ig,">=")
        arrset(pb.cnames,ig,"CON0")
        push!(irA,ig)
        push!(icA,ix_["X2"])
        push!(valA,Float64(1.5))
        push!(irA,ig)
        push!(icA,ix_["X3"])
        push!(valA,Float64(-1.0))
        for I = Int64(v_["1"]):Int64(v_["20"])
            ig,ig_,_ = s2mpj_ii("CON"*string(I),ig_)
            arrset(gtype,ig,">=")
            arrset(pb.cnames,ig,"CON"*string(I))
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
        pbm.gconst[ig_["OBJ"]] = Float64(0.1365)
        pbm.gconst[ig_["CON1"]] = Float64(213.1)
        for I = Int64(v_["2"]):Int64(v_["17"])
            pbm.gconst[ig_["CON"*string(I)]] = Float64(v_["A"*string(I)])
        end
        pbm.gconst[ig_["CON19"]] = Float64(-21.0)
        pbm.gconst[ig_["CON20"]] = Float64(110.6)
        #%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange = Vector{Float64}(undef,ngrp)
        grange[gegrps,1] = fill(Inf,pb.nge)
        for I = Int64(v_["2"]):Int64(v_["17"])
            v_["DIF"] = v_["B"*string(I)]-v_["A"*string(I)]
            arrset(grange,ig_["CON"*string(I)],Float64(v_["DIF"]))
        end
        #%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(Float64,pb.n)
        pb.xupper =    fill(Inf,pb.n)
        pb.xlower[ix_["X1"]] = 7.044148e+2
        pb.xupper[ix_["X1"]] = 9.063855e+2
        pb.xlower[ix_["X2"]] = 6.86e+1
        pb.xupper[ix_["X2"]] = 2.8888e+2
        pb.xupper[ix_["X3"]] = 1.3475e+2
        pb.xlower[ix_["X4"]] = 1.930e+2
        pb.xupper[ix_["X4"]] = 2.870966e+2
        pb.xlower[ix_["X5"]] = 2.50e+1
        pb.xupper[ix_["X5"]] = 8.41988e+1
        #%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = zeros(Float64,pb.n)
        pb.y0 = zeros(Float64,pb.m)
        if haskey(ix_,"X1")
            pb.x0[ix_["X1"]] = Float64(9.0e+2)
        else
            pb.y0[findfirst(x->x==ig_["X1"],pbm.congrps)] = Float64(9.0e+2)
        end
        if haskey(ix_,"X2")
            pb.x0[ix_["X2"]] = Float64(8.0e+1)
        else
            pb.y0[findfirst(x->x==ig_["X2"],pbm.congrps)] = Float64(8.0e+1)
        end
        if haskey(ix_,"X3")
            pb.x0[ix_["X3"]] = Float64(1.15e+2)
        else
            pb.y0[findfirst(x->x==ig_["X3"],pbm.congrps)] = Float64(1.15e+2)
        end
        if haskey(ix_,"X4")
            pb.x0[ix_["X4"]] = Float64(2.67e+2)
        else
            pb.y0[findfirst(x->x==ig_["X4"],pbm.congrps)] = Float64(2.67e+2)
        end
        if haskey(ix_,"X5")
            pb.x0[ix_["X5"]] = Float64(2.7e+1)
        else
            pb.y0[findfirst(x->x==ig_["X5"],pbm.congrps)] = Float64(2.7e+1)
        end
        #%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_  = Dict{String,Int}()
        elftv = Vector{Vector{String}}()
        it,iet_,_ = s2mpj_ii( "eY", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        elftp = Vector{Vector{String}}()
        loaset(elftp,it,1,"PI")
        it,iet_,_ = s2mpj_ii( "eC", iet_)
        loaset(elftv,it,1,"X1")
        loaset(elftv,it,2,"X2")
        loaset(elftv,it,3,"X3")
        loaset(elftv,it,4,"X4")
        loaset(elftv,it,5,"X5")
        loaset(elftp,it,1,"PI")
        #%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_      = Dict{String,Int}()
        ielftype = Vector{Int64}()
        for I = Int64(v_["1"]):Int64(v_["19"])
            v_["PI"] = Float64(I)
            v_["PI"] = 0.01+v_["PI"]
            ename = "C"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eC")
            arrset(ielftype,ie,iet_["eC"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="PI",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["PI"]))
        end
        for I = Int64(v_["1"]):Int64(v_["20"])
            v_["PI"] = Float64(I)
            v_["PI"] = 0.01+v_["PI"]
            ename = "Y"*string(I)
            ie,ie_,_  = s2mpj_ii(ename,ie_)
            arrset(pbm.elftype,ie,"eY")
            arrset(ielftype,ie,iet_["eY"])
            vname = "X1"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X1",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X2"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X2",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X3"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X3",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X4"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X4",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            vname = "X5"
            iv,ix_,pb = s2mpj_nlx(vname,ix_,pb,1,nothing,nothing,nothing)
            posev = findfirst(x->x=="X5",elftv[ielftype[ie]])
            loaset(pbm.elvar,ie,posev,iv)
            posep = findfirst(x->x=="PI",elftp[ielftype[ie]])
            loaset(pbm.elpar,ie,posep,Float64(v_["PI"]))
        end
        #%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        for ig in 1:ngrp
            arrset(pbm.grelt,ig,Int64[])
        end
        nlc = Int64[]
        ig = ig_["OBJ"]
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y17"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(-5.843e-7))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y14"])
        loaset(pbm.grelw,ig,posel,Float64(1.17e-4))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y13"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(2.358e-5))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y16"])
        loaset(pbm.grelw,ig,posel,Float64(1.502e-6))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["Y12"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(0.0321))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["Y5"])
        loaset(pbm.grelw,ig,posel,Float64(0.00423))
        posel = length(pbm.grelt[ig])+1
        loaset(pbm.grelt,ig,posel,ie_["C18"])
        arrset(nlc,length(nlc)+1,ig)
        loaset(pbm.grelw,ig,posel,Float64(1.0e-4))
        posel = posel+1
        loaset(pbm.grelt,ig,posel,ie_["C19"])
        loaset(pbm.grelw,ig,posel,Float64(37.48))
        for I = Int64(v_["1"]):Int64(v_["20"])
            ig = ig_["CON"*string(I)]
            posel = length(pbm.grelt[ig])+1
            loaset(pbm.grelt,ig,posel,ie_["Y"*string(I)])
            arrset(nlc,length(nlc)+1,ig)
            loaset(pbm.grelw,ig,posel,1.)
        end
        #%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
#    Solution
# LO SOLTN               -1.90513375
        #%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n)
        #%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        #%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower = -1*fill(Inf,pb.m)
        pb.cupper =    fill(Inf,pb.m)
        pb.clower[pb.nle+pb.neq+1:pb.m] = zeros(Float64,pb.nge)
        pb.cupper[1:pb.nge] = grange[gegrps]
        #%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons = findall(x-> x in setdiff( pbm.congrps,nlc),pbm.congrps)
        pb.pbclass = "C-COOI2-MN-5-21"
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

    elseif action == "eY"

        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        I       = round(Int,pbm.elpar[iel_][1])
        C, Y    = pbm.call("extfunc", args[1] )
        f_      = Y[1,I]
        if nargout > 1
            g_    = zeros(Float64,5)
            g_[1] = Y[2,I]
            g_[2] = Y[3,I]
            g_[3] = Y[4,I]
            g_[4] = Y[5,I]
            g_[5] = Y[6,I]
            if nargout > 2 
                H_      = zeros(Float64,5,5)
                H_[1,1] = Y[7,I]
                H_[1,2] = Y[8,I]
                H_[2,1] = H_[1,2]
                H_[2,2] = Y[9,I]
                H_[1,3] = Y[10,I]
                H_[3,1] = H_[1,3]
                H_[2,3] = Y[11,I]
                H_[3,2] = H_[2,3]
                H_[3,3] = Y[12,I]
                H_[1,4] = Y[13,I]
                H_[4,1] = H_[1,4]
                H_[2,4] = Y[14,I]
                H_[4,2] = H_[2,4]
                H_[3,4] = Y[15,I]
                H_[4,3] = H_[3,4]
                H_[4,4] = Y[16,I]
                H_[1,5] = Y[17,I]
                H_[5,1] = H_[1,5]
                H_[2,5] = Y[18,I]
                H_[5,2] = H_[2,5]
                H_[3,5] = Y[19,I]
                H_[5,3] = H_[3,5]
                H_[4,5] = Y[20,I]
                H_[5,4] = H_[4,5]
                H_[5,5] = Y[21,I]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end

    elseif action == "eC"

        EV_     = args[1]
        iel_    = args[2]
        nargout = args[3]
        pbm     = args[4]
        I       = round(Int,pbm.elpar[iel_][1])
        C, Y    = pbm.call("extfunc",args[1] )
        f_      = C[1,I]
        if nargout>1 
            g_    = zeros(Float64,5)
            g_[1] = C[2,I]
            g_[2] = C[3,I]
            g_[3] = C[4,I]
            g_[4] = C[5,I]
            g_[5] = C[6,I]
            if nargout>2 
                H_      = zeros(Float64,5,5)
                H_[1,1] = C[7,I]
                H_[1,2] = C[8,I]
                H_[2,1] = H_[1,2]
                H_[2,2] = C[9,I]
                H_[1,3] = C[10,I]
                H_[3,1] = H_[1,3]
                H_[2,3] = C[11,I]
                H_[3,2] = H_[2,3]
                H_[3,3] = C[12,I]
                H_[1,4] = C[13,I]
                H_[4,1] = H_[1,4]
                H_[2,4] = C[14,I]
                H_[4,2] = H_[2,4]
                H_[3,4] = C[15,I]
                H_[4,3] = H_[3,4]
                H_[4,4] = C[16,I]
                H_[1,5] = C[17,I]
                H_[5,1] = H_[1,5]
                H_[2,5] = C[18,I]
                H_[5,2] = H_[2,5]
                H_[3,5] = C[19,I]
                H_[5,3] = H_[3,5]
                H_[4,5] = C[20,I]
                H_[5,4] = H_[4,5]
                H_[5,5] = C[21,I]
            end
        end
        if nargout == 1
            return f_
        elseif nargout == 2
            return f_,g_
        elseif nargout == 3
            return f_,g_,H_
        end
        
    elseif action == "extfunc"

        X  = args[1]
        X1 = X[1]
        X2 = X[2]
        X3 = X[3]
        X4 = X[4]
        X5 = X[5]

        C = zeros(21,19)
        Y = zeros(21,20)

        # Function Y1
        
        A      = 4.16e+1
        Y[1,1] = X2 + X3 + A
        Y[3,1] = 1.0
        Y[4,1] = 1.0

        # Function C1
        
        A = 2.4e-2
        B = 4.62
        C[1,1] = A * X4 - B
        C[5,1] = A

        # Function Y2
        
        A       = 1.2e+1
        B       = 1.25e+1
        Y[1,2]  = A + B / C[1,1]
        Y[5,2]  = - B * C[5,1] / C[1,1]^2
        Y[16,2] = - B * C[16,1] / C[1,1]^2 + 2.0 * B * C[5,1]^2 / C[1,1]^3

        # Function C2
        
        A       = 3.535e-4
        B       = 5.311e-1
        D       = 8.705e-2
        C[1,2]  = A * X1^2 + B * X1 + D * X1 * Y[1,2]
        C[2,2]  = 2.0 * A * X1 + B + D * Y[1,2]
        C[5,2]  = D * X1 * Y[5,2]
        C[7,2]  = 2.0 * A
        C[13,2] = D * Y[5,2]
        C[16,2] = D * X1 * Y[16,2]

        # Function C3
        
        A       = 5.2e-2
        B       = 7.8e+1
        D       = 2.377e-3
        C[1,3]  = A * X1 + B + D * X1 * Y[1,2]
        C[2,3]  = A + D * Y[1,2]
        C[5,3]  = D * X1 * Y[5,2]
        C[13,3] = D * Y[5,2]
        C[16,3] = D * X1 * Y[16,2]

        # Function Y3
        
        Y[1,3]  = C[1,2] / C[1,3]
        Y[2,3]  = C[2,2] / C[1,3] - C[1,2] * C[2,3] / C[1,3]^2
        Y[5,3]  = C[5,2] / C[1,3] - C[1,2] * C[5,3] / C[1,3]^2
        Y[7,3]  = ( C[7,2] / C[1,3] - 2.0 * C[2,2] * C[2,3] / C[1,3]^2
                  - C[1,2] * C[7,3] / C[1,3]^2 + 2.0 * C[1,2] * C[2,3]^2 / C[1,3]^3 )
        Y[13,3] = ( C[13,2] / C[1,3] - C[2,2] * C[5,3] / C[1,3]^2 - C[5,2] * C[2,3] / C[1,3]^2
                  - C[1,2] * C[13,3] / C[1,3]^2 + 2.0 * C[1,2] * C[2,3] * C[5,3] / C[1,3]^3 )
        Y[16,3] = ( C[16,2] / C[1,3] - 2.0 * C[5,2] * C[5,3] / C[1,3]^2
                  - C[1,2] * C[16,3] / C[1,3]^2 + 2.0 * C[1,2] * C[5,3]^2 / C[1,3]^3 )

        # Function Y4
        
        A       = 1.9e+1
        Y[1,4]  = A * Y[1,3]
        Y[2,4]  = A * Y[2,3]
        Y[5,4]  = A * Y[5,3]
        Y[7,4]  = A * Y[7,3]
        Y[13,4] = A * Y[13,3]
        Y[16,4] = A * Y[16,3]

        # Function C4
        
        A       = 4.782e-2
        B       = 1.956e-1
        D       = 6.376e-1
        E       = 1.594
        C[1,4]  = A * (X1 - Y[1,3]) + B * ((X1 - Y[1,3])^2) / X2 + D * Y[1,4] + E * Y[1,3]
        C[2,4]  = A * (1.0 - Y[2,3]) + 2.0 * B * (X1 - Y[1,3]) * (1.0 - Y[2,3]) / X2 + D * Y[2,4] + E * Y[2,3]
        C[3,4]  = - B * ((X1 - Y[1,3])^2) / X2^2
        C[5,4]  = A * (-Y[5,3]) + 2.0 * B * (X1 - Y[1,3]) * (-Y[5,3]) / X2 + D * Y[5,4] + E * Y[5,3]
        C[7,4]  = ( A * (-Y[7,3]) + 2.0 * B * (1.0 - Y[2,3])^2 / X2
                  + 2.0 * B * (X1 - Y[1,3]) * (-Y[7,3]) / X2 + D * Y[7,4] + E * Y[7,3] )
        C[8,4]  = - 2.0 * B * (X1 - Y[1,3]) * (1.0 - Y[2,3]) / X2^2
        C[9,4]  =  2.0 * B * ((X1 - Y[1,3])^2) / X2^3
        C[13,4] = ( A * (-Y[13,3]) + 2.0 * B * (-Y[5,3]) * (1.0 - Y[2,3]) / X2
                  + 2.0 * B * (X1 - Y[1,3]) * (-Y[13,3]) / X2 + D * Y[13,4] + E * Y[13,3] )
        C[14,4] = - 2.0 * B * (X1 - Y[1,3]) * (-Y[5,3]) / X2^2
        C[16,4] = ( A * (-Y[16,3]) + 2.0 * B * (-Y[5,3])^2 / X2 + 2.0 * B * (X1 - Y[1,3]) * (-Y[16,3]) / X2
                  + D * Y[16,4] + E * Y[16,3] )

        # Function C5
        
        A      = 1.0e+2
        C[1,5] = A * X2
        C[3,5] = A

        # Function C6
        
        C[1,6]  = X1 - Y[1,3] - Y[1,4]
        C[2,6]  = 1.0 - Y[2,3] - Y[2,4]
        C[5,6]  = - Y[5,3] - Y[5,4]
        C[7,6]  = - Y[7,3] - Y[7,4]
        C[13,6] = - Y[13,3] - Y[13,4]
        C[16,6] = - Y[16,3] - Y[16,4]

        # Function C7
        
        A       = 9.5e-1
        C[1,7]  = A - C[1,4] / C[1,5]
        C[2,7]  = - C[2,4] / C[1,5]
        C[3,7]  = - C[3,4] / C[1,5] + C[1,4] * C[3,5] / C[1,5]^2
        C[5,7]  = - C[5,4] / C[1,5]
        C[7,7]  = - C[7,4] / C[1,5]
        C[8,7]  = - C[8,4] / C[1,5] + C[2,4] * C[3,5] / C[1,5]^2
        C[9,7]  = ( - C[9,4] / C[1,5] + 2.0 * C[3,4] * C[3,5] / C[1,5]^2
                  + C[1,4] * C[9,5] / C[1,5]^2 - 2.0 * C[1,4] * C[3,5]^2 / C[1,5]^3 )
        C[13,7] = - C[13,4] / C[1,5]
        C[14,7] = - C[14,4] / C[1,5] + C[5,4] * C[3,5] / C[1,5]^2 + C[1,4] * C[14,5] / C[1,5]^2
        C[16,7] = - C[16,4] / C[1,5]

        # Function Y5
        
        Y[1,5]  = C[1,6]  * C[1,7]
        Y[2,5]  = C[2,6]  * C[1,7] + C[1,6] * C[2,7]
        Y[3,5]  = C[1,6]  * C[3,7]
        Y[5,5]  = C[5,6]  * C[1,7] + C[1,6] * C[5,7]
        Y[7,5]  = C[7,6]  * C[1,7] + 2.0 * C[2,6] * C[2,7] + C[1,6] * C[7,7]
        Y[8,5]  = C[2,6]  * C[3,7] + C[1,6] * C[8,7]
        Y[9,5]  = C[1,6]  * C[9,7]
        Y[13,5] = C[13,6] * C[1,7] + C[5,6] * C[2,7] + C[2,6] * C[5,7] + C[1,6] * C[13,7]
        Y[14,5] = C[5,6]  * C[3,7] + C[1,6] * C[14,7]
        Y[16,5] = C[16,6] * C[1,7] + 2.0 * C[5,6] * C[5,7] + C[1,6] * C[16,7]

        # Function Y6
        
        Y[1,6]  = X1 - Y[1,3] - Y[1,4] - Y[1,5]
        Y[2,6]  = 1.0 - Y[2,3] - Y[2,4] - Y[2,5]
        Y[3,6]  = - Y[3,5]
        Y[5,6]  = - Y[5,3] - Y[5,4] - Y[5,5]
        Y[7,6]  = - Y[7,3] - Y[7,4] - Y[7,5]
        Y[8,6]  = - Y[8,5]
        Y[9,6]  = - Y[9,5]
        Y[13,6] = - Y[13,3] - Y[13,4] - Y[13,5]
        Y[14,6] = - Y[14,5]
        Y[16,6] = - Y[16,3] - Y[16,4] - Y[16,5]

        # Function C8
        
        A       = 9.95e-1
        C[1,8]  = A * (Y[1,4] + Y[1,5])
        C[2,8]  = A * (Y[2,4] + Y[2,5])
        C[3,8]  = A * Y[3,5]
        C[5,8]  = A * (Y[5,4] + Y[5,5])
        C[7,8]  = A * (Y[7,4] + Y[7,5])
        C[8,8]  = A * Y[8,5]
        C[9,8]  = A * Y[9,5]
        C[13,8] = A * (Y[13,4] + Y[13,5])
        C[14,8] = A * Y[14,5]
        C[16,8] = A * (Y[16,4] + Y[16,5])

        # Function Y7
        
        Y[1,7] = C[1,8] / Y[1,1]
        Y[2,7] = C[2,8] / Y[1,1]
        Y[3,7] = C[3,8] / Y[1,1] - C[1,8] * Y[3,1] / Y[1,1]^2
        Y[4,7] = - C[1,8] * Y[4,1] / Y[1,1]^2
        Y[5,7] = C[5,8] / Y[1,1]
        Y[7,7] = C[7,8] / Y[1,1]
        Y[8,7] = C[8,8] / Y[1,1] - C[2,8] * Y[3,1] / Y[1,1]^2
        Y[9,7] = C[9,8] / Y[1,1] - 2.0 * C[3,8] * Y[3,1] / Y[1,1]^2 + 2.0 * C[1,8] * Y[3,1]^2 / Y[1,1]^3
        Y[10,7] = - C[2,8] * Y[4,1] / Y[1,1]^2
        Y[11,7] = ( C[11,8] / Y[1,1] - C[3,8] * Y[4,1] / Y[1,1]^2 - C[4,8] * Y[3,1] / Y[1,1]^2
                  + 2.0 * C[1,8] * Y[3,1] * Y[4,1] / Y[1,1]^3 )
        Y[12,7] = - C[4,8] * Y[4,1] / Y[1,1]^2 + 2.0 * C[1,8] * Y[4,1]^2 / Y[1,1]^3
        Y[13,7] = C[13,8] / Y[1,1]
        Y[14,7] = C[14,8] / Y[1,1] - C[5,8] * Y[3,1] / Y[1,1]^2
        Y[15,7] = - C[5,8] * Y[4,1] / Y[1,1]^2 + 2.0 * C[1,8] * Y[4,1] * Y[5,1] / Y[1,1]^3
        Y[16,7] = C[16,8] / Y[1,1]

        # Function Y8
        
        A       = 3.798e+3
        Y[1,8]  = C[1,8]  / A
        Y[2,8]  = C[2,8]  / A
        Y[3,8]  = C[3,8]  / A
        Y[5,8]  = C[5,8]  / A
        Y[7,8]  = C[7,8]  / A
        Y[8,8]  = C[8,8]  / A
        Y[9,8]  = C[9,8]  / A
        Y[13,8] = C[13,8] / A
        Y[14,8] = C[14,8] / A
        Y[16,8] = C[16,8] / A

        # Function C9
        
        A       = 6.63e-2
        B       = 3.153e-1
        C[1,9]  = Y[1,7] - A * Y[1,7] / Y[1,8] - B
        C[2,9]  = Y[2,7] - A * Y[2,7] / Y[1,8] + A * Y[1,7] * Y[2,8] / Y[1,8]^2
        C[3,9]  = Y[3,7] - A * Y[3,7] / Y[1,8] + A * Y[1,7] * Y[3,8] / Y[1,8]^2
        C[4,9]  = Y[4,7] - A * Y[4,7] / Y[1,8] + A * Y[1,7] * Y[4,8] / Y[1,8]^2
        C[5,9]  = Y[5,7] - A * Y[5,7] / Y[1,8] + A * Y[1,7] * Y[5,8] / Y[1,8]^2
        C[7,9]  = ( Y[7,7] - A * Y[7,7] / Y[1,8] + 2.0 * A * Y[2,7] * Y[2,8] / Y[1,8]^2
                  + A * Y[1,7] * Y[7,8] / Y[1,8]^2 - 2.0 * A * Y[1,7] * Y[2,8]^2 / Y[1,8]^3 )
        C[8,9]  = ( Y[8,7] - A * Y[8,7] / Y[1,8] + A * Y[2,7] * Y[3,8] / Y[1,8]^2
                  + A * Y[3,7] * Y[2,8] / Y[1,8]^2 + A * Y[1,7] * Y[8,8] / Y[1,8]^2
                  - 2.0 * A * Y[1,7] * Y[2,8] * Y[3,8] / Y[1,8]^3 )
        C[9,9]  = ( Y[9,7] - A * Y[9,7] / Y[1,8] + 2.0 * A * Y[3,7] * Y[3,8] / Y[1,8]^2
                  + A * Y[1,7] * Y[9,8] / Y[1,8]^2 - 2.0 * A * Y[1,7] * Y[3,8]^2 / Y[1,8]^3 )
        C[10,9] = ( Y[10,7] - A * Y[10,7] / Y[1,8] + A * Y[2,7] * Y[4,8]/ Y[1,8]^2
                  + A * Y[4,7] * Y[2,8] / Y[1,8]^2 + A * Y[1,7] * Y[10,8] / Y[1,8]^2
                  - 2.0 * A * Y[1,7] * Y[2,8] * Y[4,8] / Y[1,8]^3 )
        C[11,9] = ( Y[11,7] - A * Y[11,7] / Y[1,8] + A * Y[3,7] * Y[4,8]/ Y[1,8]^2
                  + A * Y[4,7] * Y[3,8] / Y[1,8]^2 + A * Y[1,7] * Y[11,8] / Y[1,8]^2 -
                  2.0 * A * Y[1,7] * Y[3,8] * Y[4,8] / Y[1,8]^3 )
        C[12,9] = ( Y[12,7] - A * Y[12,7] / Y[1,8] + 2.0 * A * Y[4,7] * Y[4,8] / Y[1,8]^2
                  + A * Y[1,7] * Y[12,8] / Y[1,8]^2 - 2.0 * A * Y[1,7] * Y[4,8]^2 / Y[1,8]^3 )
        C[13,9] = ( Y[13,7] - A * Y[13,7] / Y[1,8] + A * Y[2,7] * Y[5,8]/ Y[1,8]^2
                  + A * Y[5,7] * Y[2,8] / Y[1,8]^2 + A * Y[1,7] * Y[13,8] / Y[1,8]^2
                  - 2.0 * A * Y[1,7] * Y[2,8] * Y[5,8] / Y[1,8]^3 )
        C[14,9] = ( Y[14,7] - A * Y[14,7] / Y[1,8] + A * Y[3,7] * Y[5,8]/ Y[1,8]^2
                  + A * Y[5,7] * Y[3,8] / Y[1,8]^2 + A * Y[1,7] * Y[14,8] / Y[1,8]^2
                  - 2.0 * A * Y[1,7] * Y[3,8] * Y[5,8] / Y[1,8]^3 )
        C[15,9] = ( Y[15,7] - A * Y[15,7] / Y[1,8] + A * Y[4,7] * Y[5,8]/ Y[1,8]^2
                  + A * Y[5,7] * Y[4,8] / Y[1,8]^2 + A * Y[1,7] * Y[15,8] / Y[1,8]^2
                  - 2.0 * A * Y[1,7] * Y[4,8] * Y[5,8] / Y[1,8]^3 )
        C[16,9] = ( Y[16,7] - A * Y[16,7] / Y[1,8] + 2.0 * A * Y[5,7] * Y[5,8] / Y[1,8]^2
                  + A * Y[1,7] * Y[16,8] / Y[1,8]^2 - 2.0 * A * Y[1,7] * Y[5,8]^2 / Y[1,8]^3 )

        # Function Y9
        
        A       = 9.682e+1
        B       = 3.21e-1
        Y[1,9]  = A / C[1,9] + B * Y[1,1]
        Y[2,9]  = - A * C[2,9] / C[1,9]^2
        Y[3,9]  = - A * C[3,9] / C[1,9]^2  + B * Y[3,1]
        Y[4,9]  = - A * C[4,9] / C[1,9]^2  + B * Y[4,1]
        Y[5,9]  = - A * C[5,9] / C[1,9]^2
        Y[7,9]  = - A * C[7,9] / C[1,9]^2  + 2.0 * A * C[2,9]^2 / C[1,9]^3
        Y[8,9]  = - A * C[8,9] / C[1,9]^2  + 2.0 * A * C[2,9] * C[3,9] / C[1,9]^3
        Y[9,9]  = - A * C[9,9] / C[1,9]^2  + 2.0 * A * C[3,9]^2 / C[1,9]^3
        Y[10,9] = - A * C[10,9] / C[1,9]^2 + 2.0 * A * C[2,9] * C[4,9] / C[1,9]^3
        Y[11,9] = - A * C[11,9] / C[1,9]^2 + 2.0 * A * C[3,9] * C[4,9] / C[1,9]^3
        Y[12,9] = - A * C[12,9] / C[1,9]^2 + 2.0 * A * C[4,9]^2 / C[1,9]^3
        Y[13,9] = - A * C[13,9] / C[1,9]^2 + 2.0 * A * C[2,9] * C[5,9] / C[1,9]^3
        Y[14,9] = - A * C[14,9] / C[1,9]^2 + 2.0 * A * C[3,9] * C[5,9] / C[1,9]^3
        Y[15,9] = - A * C[15,9] / C[1,9]^2 + 2.0 * A * C[4,9] * C[5,9] / C[1,9]^3
        Y[16,9] = - A * C[16,9] / C[1,9]^2 + 2.0 * A * C[5,9]^2 / C[1,9]^3

        # Function Y10
        
        A        = 2.29
        B        = 1.258
        D        = 1.29
        E        = 1.71
        Y[1,10]  = A * Y[1,3]  + B * Y[1,4]   + D * Y[1,5] + E * Y[1,6]
        Y[2,10]  = A * Y[2,3]  + B * Y[2,4]   + D * Y[2,5] + E * Y[2,6]
        Y[3,10]  = D * Y[3,5]  + E * Y[3,6]
        Y[5,10]  = A * Y[5,3]  + B * Y[5,4]   + D * Y[5,5] + E * Y[5,6]
        Y[7,10]  = A * Y[7,3]  + B * Y[7,4]   + D * Y[7,5] + E * Y[7,6]
        Y[8,10]  = D * Y[8,5]  + E * Y[8,6]
        Y[9,10]  = D * Y[9,5]  + E * Y[9,6]
        Y[13,10] = A * Y[13,3] + B * Y[13,4] + D * Y[13,5] + E * Y[13,6]
        Y[14,10] = D * Y[14,5] + E * Y[14,6]
        Y[16,10] = A * Y[16,3] + B * Y[16,4] + D * Y[16,5] + E * Y[16,6]

        # Function Y11
        
        A        = 1.71
        B        = 5.8e-1
        D        = 4.52e-1
        Y[1,11]  = A * X1 + B * Y[1,3] - D * Y[1,4]
        Y[2,11]  = A + B * Y[2,3] - D * Y[2,4]
        Y[5,11]  = B * Y[5,3]  - D * Y[5,4]
        Y[7,11]  = B * Y[7,3]  - D * Y[7,4]
        Y[13,11] = B * Y[13,3] - D * Y[13,4]
        Y[16,11] = B * Y[16,3] - D * Y[16,4]

        # Function C10
        
        A       = 1.23e+1
        B       = 7.523e+2
        C[1,10] = A / B

        # Function C11
        
        A        = 1.74125
        C[1,11]  = A * X1 * Y[1,2]
        C[2,11]  = A * Y[1,2]
        C[5,11]  = A * X1 * Y[5,2]
        C[13,11] = A * Y[5,2]
        C[16,11] = A * X1 * Y[16,2]

        # Function C12
        
        A        = 9.995e-1
        B        = 1.998e+3
        C[1,12]  = A * Y[1,10] + B
        C[2,12]  = A * Y[2,10]
        C[3,12]  = A * Y[3,10]
        C[5,12]  = A * Y[5,10]
        C[7,12]  = A * Y[7,10]
        C[8,12]  = A * Y[8,10]
        C[9,12]  = A * Y[9,10]
        C[13,12] = A * Y[13,10]
        C[14,12] = A * Y[14,10]
        C[16,12] = A * Y[16,10]

        # Function Y12
        
        Y[1,12]  = C[1,10] * X1 + C[1,11] / C[1,12]
        Y[2,12]  = C[1,10] + C[2,11] / C[1,12] - C[1,11] * C[2,12] / C[1,12]^2
        Y[3,12]  = - C[1,11] * C[3,12] / C[1,12]^2
        Y[5,12]  = C[5,11] / C[1,12] - C[1,11] * C[5,12] / C[1,12]^2
        Y[7,12]  = ( C[7,11] / C[1,12] - 2.0 * C[2,11] * C[2,12] / C[1,12]^2
                   - C[1,11] * C[7,12] / C[1,12]^2 + 2.0 * C[1,11] * C[2,12]^2 / C[1,12]^3 )
        Y[8,12]  = ( C[2,11] * C[3,12] / C[1,12]^2 - C[1,11] * C[8,12] / C[1,12]^2
                   + 2.0 * C[1,11] * C[2,12] * C[3,12] / C[1,12]^3 )
        Y[9,12]  = - C[1,11] * C[9,12] / C[1,12]^2 + 2.0 * C[1,11] * C[3,12]^2 / C[1,12]^3
        Y[13,12] = ( C[13,11] / C[1,12] - C[2,11] * C[5,12] / C[1,12]^2 - C[5,11] * C[2,12] / C[1,12]^2
                   - C[1,11] * C[13,12] / C[1,12]^2 + 2.0 * C[1,11] * C[2,12] * C[5,12] / C[1,12]^3 )
        Y[14,12] = ( - C[5,11] * C[3,12] / C[1,12]^2 - C[1,11] * C[14,12] / C[1,12]^2
                   + 2.0 * C[1,11] * C[3,12] * C[5,12] / C[1,12]^3 )
        Y[16,12] = ( C[16,11] / C[1,12] - 2.0 * C[5,11] * C[5,12] / C[1,12]^2
                   - C[1,11] * C[16,12] / C[1,12]^2 + 2.0 * C[1,11] * C[5,12]^2 / C[1,12]^3 )

        # Function Y13
        
        A        = 1.75
        Y[1,13]  = C[1,12] - A * Y[1,2]
        Y[2,13]  = C[2,12]
        Y[3,13]  = C[3,12]
        Y[5,13]  = C[5,12] - A * Y[5,2]
        Y[7,13]  = C[7,12]
        Y[8,13]  = C[8,12]
        Y[9,13]  = C[9,12]
        Y[13,13] = C[13,12]
        Y[14,13] = C[14,12]
        Y[16,13] = C[16,12] - A * Y[16,2]

        # Function Y14
        
        A        = 3.623e+3
        B        = 6.44e+1
        D        = 1.46312e+5
        F        = 5.84e+1
        Y[1,14]  = A + B * X2 + F * X3 + D / ( Y[1,9] + X5)
        Y[2,14]  = - D * Y[2,9] / ( Y[1,9] + X5)^2
        Y[3,14]  = B - D * Y[3,9] / ( Y[1,9] + X5)^2
        Y[4,14]  = F - D * Y[4,9] / ( Y[1,9] + X5)^2
        Y[5,14]  = - D * Y[5,9] / ( Y[1,9] + X5)^2
        Y[6,14]  = - D / ( Y[1,9] + X5)^2
        Y[7,14]  = - D * Y[7,9] / ( Y[1,9] + X5)^2 + 2.0 * D * Y[2,9]^2 / (Y[1,9] + X5)^3
        Y[8,14]  = - D * Y[8,9] / ( Y[1,9] + X5)^2 + 2.0 * D * Y[2,9] * Y[3,9] / (Y[1,9] + X5)^3
        Y[9,14]  = - D * Y[9,9] / ( Y[1,9] + X5)^2 + 2.0 * D * Y[3,9]^2 / (Y[1,9] + X5)^3
        Y[10,14] = - D * Y[10,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[2,9]*Y[4,9]/(Y[1,9] + X5)^3
        Y[11,14] = - D * Y[11,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[3,9]*Y[4,9]/(Y[1,9] + X5)^3
        Y[12,14] = - D * Y[12,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[4,9]^2 / (Y[1,9] + X5)^3
        Y[13,14] = - D * Y[13,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[2,9]*Y[5,9]/(Y[1,9] + X5)^3
        Y[14,14] = - D * Y[14,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[3,9]*Y[5,9]/(Y[1,9] + X5)^3
        Y[15,14] = - D * Y[15,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[4,9]*Y[5,9]/(Y[1,9] + X5)^3
        Y[16,14] = - D * Y[16,9]/(Y[1,9]+X5)^2 + 2.0 * D * Y[5,9]^2/(Y[1,9] + X5)^3
        Y[17,14] = 2.0 * D * Y[2,9]/(Y[1,9] + X5)^3
        Y[18,14] = 2.0 * D * Y[3,9]/(Y[1,9] + X5)^3
        Y[19,14] = 2.0 * D * Y[4,9]/(Y[1,9] + X5)^3
        Y[20,14] = 2.0 * D * Y[5,9]/(Y[1,9] + X5)^3
        Y[21,14] = 2.0 * D / (Y[1,9] + X5)^3

        # Function C13
        
        A        = 9.995e-1
        B        = 1.121e-1
        D        = 4.8e+1
        E        = 5.095e+3
        F        = 6.08e+1
        C[1,13]  = A * Y[1,10] - B * Y[1,14] + F * X2 + D * X4 - E
        C[2,13]  = A * Y[2,10] - B * Y[2,14]
        C[3,13]  = F + A * Y[3,10] - B * Y[3,14]
        C[4,13]  = - B * Y[4,14]
        C[5,13]  = D + A * Y[5,10] - B * Y[5,14]
        C[6,13]  = - B * Y[6,14]
        C[7,13]  = A * Y[7,10] - B * Y[7,14]
        C[8,13]  = A * Y[8,10] - B * Y[8,14]
        C[9,13]  = A * Y[9,10] - B * Y[9,14]
        C[10,13] = - B * Y[10,14]
        C[11,13] = - B * Y[11,14]
        C[12,13] = - B * Y[12,14]
        C[13,13] = A * Y[13,10] - B * Y[13,14]
        C[14,13] = A * Y[14,10] - B * Y[14,14]
        C[15,13] = - B * Y[15,14]
        C[16,13] = A * Y[16,10] - B * Y[16,14]
        C[17,13] = - B * Y[17,14]
        C[18,13] = - B * Y[18,14]
        C[19,13] = - B * Y[19,14]
        C[20,13] = - B * Y[20,14]
        C[21,13] = - B * Y[21,14]

        # Function Y15
        
        Y[1,15]  = Y[1,13] / C[1,13]
        Y[2,15]  = Y[2,13] / C[1,13] - Y[1,13] * C[2,13] / C[1,13]^2
        Y[3,15]  = Y[3,13] / C[1,13] - Y[1,13] * C[3,13] / C[1,13]^2
        Y[4,15]  = - Y[1,13] * C[4,13] / C[1,13]^2
        Y[5,15]  = Y[5,13] / C[1,13] - Y[1,13] * C[5,13] / C[1,13]^2
        Y[6,15]  = - Y[1,13] * C[6,13] / C[1,13]^2
        Y[7,15]  = ( Y[7,13] / C[1,13] - 2.0 * Y[2,13] * C[2,13] / C[1,13]^2 - Y[1,13] * C[7,13] / C[1,13]^2
                   + 2.0 * Y[1,13] * C[2,13]^2 / C[1,13]^3 )
        Y[8,15]  = ( Y[8,13] / C[1,13] - Y[2,13] * C[3,13] / C[1,13]^2 - Y[3,13] * C[2,13] / C[1,13]^2
                   - Y[1,13] * C[8,13] / C[1,13]^2 + 2.0 * Y[1,13] * C[2,13]*C[3,13] / C[1,13]^3 )
        Y[9,15]  = ( Y[9,13] / C[1,13] - 2.0 * Y[3,13]*C[3,13]/C[1,13]^2
                   - Y[1,13]*C[9,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[3,13]^2 / C[1,13]^3 )
        Y[10,15] = ( Y[10,13]/C[1,13] - Y[2,13]*C[4,13]/C[1,13]^2 - Y[4,13]*C[2,13]/C[1,13]^2
                   - Y[1,13]*C[10,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[2,13]*C[4,13]/C[1,13]^3 )
        Y[11,15] = ( Y[3,13]*C[4,13]/C[1,13]^2 - Y[1,13]*C[11,13]/C[1,13]^2
                   + 2.0 * Y[1,13]*C[3,13]*C[4,13]/C[1,13]^3 )
        Y[12,15] = - Y[1,13]*C[12,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[4,13]^2 / C[1,13]^3
        Y[13,15] = ( Y[13,13]/C[1,13] - Y[2,13]*C[5,13]/C[1,13]^2 - Y[5,13]*C[2,13]/C[1,13]^2
                   - Y[1,13]*C[13,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[2,13]*C[5,13]/C[1,13]^3 )
        Y[14,15] = ( Y[14,13]/C[1,13] - Y[3,13]*C[5,13]/C[1,13]^2 - Y[5,13]*C[3,13]/C[1,13]^2
                   - Y[1,13]*C[14,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[3,13]*C[5,13]/C[1,13]^3 )
        Y[15,15] = ( - Y[5,13]*C[4,13]/C[1,13]^2 - Y[1,13]*C[15,13]/C[1,13]^2
                   + 2.0 * Y[1,13]*C[4,13]*C[5,13]/C[1,13]^3 )
        Y[16,15] = ( Y[16,13]/C[1,13] - 2.0 * Y[5,13]*C[5,13]/C[1,13]^2
                   - Y[1,13]*C[16,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[5,13]^2 / C[1,13]^3 )
        Y[17,15] = ( - Y[2,13]*C[6,13]/C[1,13]^2 - Y[1,13]*C[17,13]/C[1,13]^2
                   + 2.0 * Y[1,13]*C[2,13]*C[6,13]/C[1,13]^3 )
        Y[18,15] = ( - Y[3,13]*C[6,13]/C[1,13]^2 - Y[1,13]*C[18,13]/C[1,13]^2
                   + 2.0 * Y[1,13]*C[3,13]*C[6,13]/C[1,13]^3 )
        Y[19,15] = - Y[1,13]*C[19,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[4,13]*C[6,13]/C[1,13]^3
        Y[20,15] = ( - Y[5,13]*C[6,13]/C[1,13]^2 - Y[1,13]*C[20,13]/C[1,13]^2
                   + 2.0 * Y[1,13]*C[5,13]*C[6,13]/C[1,13]^3 )
        Y[21,15] = - Y[1,13]*C[21,13]/C[1,13]^2 + 2.0 * Y[1,13]*C[6,13]^2 / C[1,13]^3

        # Function Y16
        
        A        = 1.48e+5
        B        = -3.31e+5
        D        = 4.0e+1
        E        = -6.1e+1
        Y[1,16]  = A + B*Y[1,15] + D*Y[1,13] + E*Y[1,15]*Y[1,13]
        Y[2,16]  = B*Y[2,15] + D*Y[2,13] + E*(Y[2,15]*Y[1,13] + Y[1,15]*Y[2,13])
        Y[3,16]  = B*Y[3,15] + D*Y[3,13] + E*(Y[3,15]*Y[1,13] + Y[1,15]*Y[3,13])
        Y[4,16]  = B*Y[4,15] + E*Y[4,15]*Y[1,13]
        Y[5,16]  = B*Y[5,15] + D*Y[5,13] + E*(Y[5,15]*Y[1,13] + Y[1,15]*Y[5,13])
        Y[6,16]  = B*Y[6,15] + E*Y[6,15]*Y[1,13]
        Y[7,16]  = B*Y[7,15] + D*Y[7,13] + E*(Y[7,15]*Y[1,13] + Y[1,15]*Y[7,13] + 2.0*Y[2,15]*Y[2,13])
        Y[8,16]  = B*Y[8,15] + D*Y[8,13] + E*(Y[8,15]*Y[1,13] + Y[1,15]*Y[8,13] + Y[2,15]*Y[3,13] + Y[3,15]*Y[2,13])
        Y[9,16]  = B*Y[9,15] + D*Y[9,13] + E*(Y[9,15]*Y[1,13] + Y[1,15]*Y[9,13] + 2.0*Y[3,15]*Y[3,13])
        Y[10,16] = B*Y[10,15] + E*(Y[10,15]*Y[1,13] + Y[4,15]*Y[2,13])
        Y[11,16] = B*Y[11,15] + E*(Y[11,15]*Y[1,13] + Y[4,15]*Y[3,13])
        Y[12,16] = B*Y[12,15] + E*Y[12,15]*Y[1,13]
        Y[13,16] = ( B*Y[13,15] + D*Y[13,13] + E*(Y[13,15]*Y[1,13] + Y[1,15]*Y[13,13]
                   + Y[2,15]*Y[5,13] + Y[5,15]*Y[2,13]) )
        Y[14,16] = ( B*Y[14,15] + D*Y[14,13] + E*(Y[14,15]*Y[1,13] + Y[1,15]*Y[14,13]
                   + Y[3,15]*Y[5,13] + Y[5,15]*Y[3,13]) )
        Y[15,16] = B*Y[15,15] + E*(Y[15,15]*Y[1,13] + Y[4,15]*Y[5,13])
        Y[16,16] = B*Y[16,15] + D*Y[16,13] + E*(Y[16,15]*Y[1,13] + Y[1,15]*Y[16,13] + 2.0*Y[5,15]*Y[5,13])
        Y[17,16] = B*Y[17,15] + E*(Y[17,15]*Y[1,13] + Y[6,15]*Y[2,13])
        Y[18,16] = B*Y[18,15] + E*(Y[18,15]*Y[1,13] + Y[6,15]*Y[3,13])
        Y[19,16] = B*Y[19,15] + E*Y[19,15]*Y[1,13]
        Y[20,16] = B*Y[20,15] + E*(Y[20,15]*Y[1,13] + Y[6,15]*Y[5,13])
        Y[21,16] = B*Y[21,15] + E*Y[21,15]*Y[1,13]

        # Function C14
        
        A       = 2.324e+3
        B       = 2.874e+7
        C[1,14] = A * Y[1,10] - B * Y[1,2]
        C[2,14] = A * Y[2,10]
        C[3,14] = A * Y[3,10]
        C[5,14] = A * Y[5,10] - B * Y[5,2]
        C[7,14] = A * Y[7,10]
        C[8,14] = A * Y[8,10]
        C[9,14] = A * Y[9,10]
        C[13,14]= A * Y[13,10]
        C[14,14]= A * Y[14,10]
        C[16,14]= A * Y[16,10] - B * Y[16,2]

        # Function Y17
        
        A        = 1.413e+7
        B        = -1.328e+3
        D        = -5.31e+2
        Y[1,17]  = A + B*Y[1,10] + D*Y[1,11] + C[1,14]/C[1,12]
        Y[2,17]  = B*Y[2,10] + D*Y[2,11] + C[2,14]/C[1,12] - C[1,14]*C[2,12]/C[1,12]^2
        Y[3,17]  = B*Y[3,10] + C[3,14]/C[1,12] - C[1,14]*C[3,12]/C[1,12]^2
        Y[5,17]  = B*Y[5,10] + D*Y[5,11] + C[5,14]/C[1,12] - C[1,14]*C[5,12]/C[1,12]^2
        Y[7,17]  = ( B*Y[7,10] + D*Y[7,11] + C[7,14]/C[1,12] - C[1,14]*C[7,12]/C[1,12]^2
                   - 2.0*C[2,14]*C[2,12]/C[1,12]^2 + 2.0*C[1,14]*C[2,12]^2/C[1,12]^3 )
        Y[8,17]  = ( B*Y[8,10] + C[8,14]/C[1,12] - C[2,14]*C[3,12]/C[1,12]^2 - C[3,14]*C[2,12]/C[1,12]^2
                   - C[1,14]*C[8,12]/C[1,12]^2 + 2.0*C[1,14]*C[2,12]*C[3,12]/C[1,12]^3 )
        Y[9,17]  = ( B*Y[9,10] + C[9,14]/C[1,12] - C[1,14]*C[9,12]/C[1,12]^2
                   - 2.0*C[3,14]*C[3,12]/C[1,12]^2 + 2.0*C[1,14]*C[3,12]^2/C[1,12]^3 )
        Y[13,17] = ( B*Y[13,10]+ D*Y[13,11]+C[13,14]/C[1,12]-C[2,14]*C[5,12]/C[1,12]^2
                   - C[5,14]*C[2,12]/C[1,12]^2 - C[1,14]*C[13,12]/C[1,12]^2 + 2.0*C[1,14]*C[2,12]*C[5,12]/C[1,12]^3 )
        Y[14,17] = ( B*Y[14,10]+C[14,14]/C[1,12]-C[3,14]*C[5,12]/C[1,12]^2 - C[5,14]*C[3,12]/C[1,12]^2
                   - C[1,14]*C[14,12]/C[1,12]^2 + 2.0*C[1,14]*C[3,12]*C[5,12]/C[1,12]^3 )
        Y[16,17] = ( B*Y[16,10]+D*Y[16,11]+C[16,14]/C[1,12]-C[1,14]*C[16,12]/C[1,12]^2
                   - 2.0*C[5,14]*C[5,12]/C[1,12]^2 + 2.0*C[1,14]*C[5,12]^2 /C[1,12]^3 )

        # Function C15
        
        A        = 5.2e-1
        C[1,15]  = Y[1,13]/Y[1,15] - Y[1,13]/A
        C[2,15]  = Y[2,13]/Y[1,15] - Y[1,13]*Y[2,15]/Y[1,15]^2 - Y[2,13]/A
        C[3,15]  = Y[3,13]/Y[1,15] - Y[1,13]*Y[3,15]/Y[1,15]^2 - Y[3,13]/A
        C[4,15]  = - Y[1,13]*Y[4,15]/Y[1,15]^2
        C[5,15]  = Y[5,13]/Y[1,15] - Y[1,13]*Y[5,15]/Y[1,15]^2 - Y[5,13]/A
        C[6,15]  = - Y[1,13]*Y[6,15]/Y[1,15]^2
        C[7,15]  = ( Y[7,13]/Y[1,15] - 2.0*Y[2,13]*Y[2,15]/Y[1,15]^2
                   - Y[1,13]*Y[7,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[2,15]^2/Y[1,15]^3 - Y[7,13]/A )
        C[8,15]  = ( Y[8,13]/Y[1,15] - Y[2,13]*Y[3,15]/Y[1,15]^2 - Y[3,13]*Y[2,15]/Y[1,15]^2
                   - Y[1,13]*Y[8,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[2,15]*Y[3,15]/Y[1,15]^3 - Y[8,13]/A )
        C[9,15]  = ( Y[9,13]/Y[1,15] - 2.0*Y[3,13]*Y[3,15]/Y[1,15]^2
                   - Y[1,13]*Y[9,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[3,15]^2/Y[1,15]^3 - Y[9,13]/A )
        C[10,15] = ( Y[10,13]/Y[1,15] - Y[2,13]*Y[4,15]/Y[1,15]^2 - Y[4,13]*Y[2,15]/Y[1,15]^2
                   - Y[1,13]*Y[10,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[2,15]*Y[4,15]/Y[1,15]^3 )
        C[11,15] = Y[3,13]*Y[4,15]/Y[1,15]^2 - Y[1,13]*Y[11,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[3,15]*Y[4,15]/Y[1,15]^3
        C[12,15] = - Y[1,13]*Y[12,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[4,15]^2/Y[1,15]^3
        C[13,15] = ( Y[13,13]/Y[1,15] - Y[2,13]*Y[5,15]/Y[1,15]^2 - Y[5,13]*Y[2,15]/Y[1,15]^2
                   - Y[1,13]*Y[13,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[2,15]*Y[5,15]/Y[1,15]^3 - Y[13,13]/A )
        C[14,15] = ( Y[14,13]/Y[1,15] - Y[3,13]*Y[5,15]/Y[1,15]^2 - Y[5,13]*Y[3,15]/Y[1,15]^2
                   - Y[1,13]*Y[14,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[3,15]*Y[5,15]/Y[1,15]^3 - Y[14,13]/A )
        C[15,15] = - Y[5,13]*Y[4,15]/Y[1,15]^2 - Y[1,13]*Y[15,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[4,15]*Y[5,15]/Y[1,15]^3
        C[16,15] = ( Y[16,13]/Y[1,15] - 2.0*Y[5,13]*Y[5,15]/Y[1,15]^2
                   - Y[1,13]*Y[16,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[5,15]^2 / Y[1,15]^3 - Y[16,13]/A )
        C[17,15] = - Y[2,13]*Y[6,15]/Y[1,15]^2 - Y[1,13]*Y[17,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[2,15]*Y[6,15]/Y[1,15]^3
        C[18,15] = - Y[3,13]*Y[6,15]/Y[1,15]^2 - Y[1,13]*Y[18,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[3,15]*Y[6,15]/Y[1,15]^3
        C[19,15] = - Y[1,13]*Y[19,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[4,15]*Y[6,15]/Y[1,15]^3
        C[20,15] = - Y[5,13]*Y[6,15]/Y[1,15]^2 - Y[1,13]*Y[20,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[5,15]*Y[6,15]/Y[1,15]^3
        C[21,15] = - Y[1,13]*Y[21,15]/Y[1,15]^2 + 2.0*Y[1,13]*Y[6,15]^2 / Y[1,15]^3

        # Function C16
        
        A       = 1.104
        B       = 7.2e-1
        C[1,16] = A - B*Y[1,15]
        for i in 2:21
            C[i,16] = - B * Y[i,15]
        end

        # Function C17
        
        C[1,17]  = Y[1,9] + X5
        C[2,17]  = Y[2,9]
        C[3,17]  = Y[3,9]
        C[4,17]  = Y[4,9]
        C[5,17]  = Y[5,9]
        C[6,17]  = 1.0
        C[7,17]  = Y[7,9]
        C[8,17]  = Y[8,9]
        C[9,17]  = Y[9,9]
        C[10,17] = Y[10,9]
        C[11,17] = Y[11,9]
        C[12,17] = Y[12,9]
        C[13,17] = Y[13,9]
        C[14,17] = Y[14,9]
        C[15,17] = Y[15,9]
        C[16,17] = Y[16,9]

        # Function C18
        
        C[1,18] = C[1,15]/C[1,16]
        C[2,18] = C[2,15]/C[1,16] - C[1,15]*C[2,16]/C[1,16]^2
        C[3,18] = C[3,15]/C[1,16] - C[1,15]*C[3,16]/C[1,16]^2
        C[4,18] = C[4,15]/C[1,16] - C[1,15]*C[4,16]/C[1,16]^2
        C[5,18] = C[5,15]/C[1,16] - C[1,15]*C[5,16]/C[1,16]^2
        C[6,18] = C[6,15]/C[1,16] - C[1,15]*C[6,16]/C[1,16]^2
        C[7,18] = ( C[7,15]/C[1,16] - 2.0*C[2,15]*C[2,16]/C[1,16]^2
                  - C[1,15]*C[7,16]/C[1,16]^2 + 2.0*C[1,15]*C[2,16]^2/C[1,16]^3 )
        C[8,18] = ( C[8,15]/C[1,16] - C[2,15]*C[3,16]/C[1,16]^2 - C[3,15]*C[2,16]/C[1,16]^2
                  - C[1,15]*C[8,16]/C[1,16]^2 + 2.0*C[1,15]*C[2,16]*C[3,16]/C[1,16]^3 )
        C[9,18] = ( C[9,15]/C[1,16] - 2.0*C[3,15]*C[3,16]/C[1,16]^2
                  - C[1,15]*C[9,16]/C[1,16]^2 + 2.0*C[1,15]*C[3,16]^2 / C[1,16]^3 )
        C[10,18]= ( C[10,15]/C[1,16] - C[2,15]*C[4,16]/C[1,16]^2 - C[4,15]*C[2,16]/C[1,16]^2
                  - C[1,15]*C[10,16]/C[1,16]^2 + 2.0*C[1,15]*C[2,16]*C[4,16]/C[1,16]^3 )
        C[11,18]= ( C[11,15]/C[1,16] - C[3,15]*C[4,16]/C[1,16]^2 - C[4,15]*C[3,16]/C[1,16]^2
                  - C[1,15]*C[11,16]/C[1,16]^2 + 2.0*C[1,15]*C[3,16]*C[4,16]/C[1,16]^3 )
        C[12,18]= ( C[12,15]/C[1,16] - 2.0*C[4,15]*C[4,16]/C[1,16]^2
                  - C[1,15]*C[12,16]/C[1,16]^2 + 2.0*C[1,15]*C[4,16]^2/C[1,16]^3 )
        C[13,18]= ( C[13,15]/C[1,16] - C[2,15]*C[5,16]/C[1,16]^2 - C[5,15]*C[2,16]/C[1,16]^2
                  - C[1,15]*C[13,16]/C[1,16]^2 + 2.0*C[1,15]*C[2,16]*C[5,16]/C[1,16]^3 )
        C[14,18]= ( C[14,15]/C[1,16] - C[3,15]*C[5,16]/C[1,16]^2 - C[5,15]*C[3,16]/C[1,16]^2
                  - C[1,15]*C[14,16]/C[1,16]^2 + 2.0*C[1,15]*C[3,16]*C[5,16]/C[1,16]^3 )
        C[15,18]= ( C[15,15]/C[1,16] - C[4,15]*C[5,16]/C[1,16]^2 - C[5,15]*C[4,16]/C[1,16]^2
                  - C[1,15]*C[15,16]/C[1,16]^2 + 2.0*C[1,15]*C[4,16]*C[5,16]/C[1,16]^3 )
        C[16,18]= ( C[16,15]/C[1,16] - 2.0*C[5,15]*C[5,16]/C[1,16]^2
                  - C[1,15]*C[16,16]/C[1,16]^2 + 2.0*C[1,15]*C[5,16]^2/C[1,16]^3 )
        C[17,18]= ( C[17,15]/C[1,16] - C[2,15]*C[6,16]/C[1,16]^2 - C[6,15]*C[2,16]/C[1,16]^2
                  - C[1,15]*C[17,16]/C[1,16]^2 + 2.0*C[1,15]*C[2,16]*C[6,16]/C[1,16]^3 )
        C[18,18]= ( C[18,15]/C[1,16] - C[3,15]*C[6,16]/C[1,16]^2 - C[6,15]*C[3,16]/C[1,16]^2
                  - C[1,15]*C[18,16]/C[1,16]^2 + 2.0*C[1,15]*C[3,16]*C[6,16]/C[1,16]^3 )
        C[19,18]= ( C[19,15]/C[1,16] - C[4,15]*C[6,16]/C[1,16]^2 - C[6,15]*C[4,16]/C[1,16]^2
                  - C[1,15]*C[19,16]/C[1,16]^2 + 2.0*C[1,15]*C[4,16]*C[6,16]/C[1,16]^3 )
        C[20,18]= ( C[20,15]/C[1,16] - C[5,15]*C[6,16]/C[1,16]^2 - C[6,15]*C[5,16]/C[1,16]^2
                  - C[1,15]*C[20,16]/C[1,16]^2 + 2.0*C[1,15]*C[5,16]*C[6,16]/C[1,16]^3 )
        C[21,18]= ( C[21,15]/C[1,16] - 2.0*C[6,15]*C[6,16]/C[1,16]^2
                  - C[1,15]*C[21,16]/C[1,16]^2 + 2.0*C[1,15]*C[6,16]^2 / C[1,16]^3 )

        # Function C19
        
        C[1,19] = Y[1,2]/C[1,12]
        C[2,19] = - Y[1,2]*C[2,12]/C[1,12]^2
        C[3,19] = - Y[1,2]*C[3,12]/C[1,12]^2
        C[5,19] = Y[5,2]/C[1,12] - Y[1,2]*C[5,12]/C[1,12]^2
        C[7,19] = Y[1,2]*(2.0*C[2,12]^2/C[1,12]^3 - C[7,12]/C[1,12]^2)
        C[8,19] = Y[1,2]*(2.0*C[2,12]*C[3,12]/C[1,12]^3 - C[8,12]/C[1,12]^2)
        C[9,19] = Y[1,2]*(2.0*C[3,12]^2/C[1,12]^3 - C[9,12]/C[1,12]^2)
        C[13,19]= Y[1,2]*(2.0*C[2,12]*C[5,12]/C[1,12]^3 - C[13,12]/C[1,12]^2) - Y[5,2]*C[2,12]/C[1,12]^2
        C[14,19]= Y[1,2]*(2.0*C[3,12]*C[5,12]/C[1,12]^3 - C[14,12]/C[1,12]^2) - Y[5,2]*C[3,12]/C[1,12]^2
        C[16,19]= ( Y[16,2]/C[1,12] - 2.0*Y[5,2]*C[5,12]/C[1,12]^2
                  - Y[1,2]*C[16,12]/C[1,12]^2 + 2.0*Y[1,2]*C[5,12]^2/C[1,12]^3 )

        # Function Y18
        
        A        = -2.8e-1 / 7.2e-1
        Y[1,18]  = Y[1,4] + A / Y[1,5]
        Y[2,18]  = Y[2,4] - A*Y[2,5]/Y[1,5]^2
        Y[3,18]  = - A*Y[3,5]/Y[1,5]^2
        Y[5,18]  = Y[5,4] - A*Y[5,5]/Y[1,5]^2
        Y[7,18]  = Y[7,4] - A*Y[7,5]/Y[1,5]^2 + 2.0*A*Y[2,5]^2 / Y[1,5]^3
        Y[8,18]  = - A*Y[8,5]/Y[1,5]^2 + 2.0*A*Y[2,5]*Y[3,5]/Y[1,5]^3
        Y[9,18]  = - A*Y[9,5]/Y[1,5]^2 + 2.0*A*Y[3,5]^2 / Y[1,5]^3
        Y[13,18] = Y[13,4]-A*Y[13,5]/Y[1,5]^2+ 2.0*A*Y[2,5]*Y[5,5]/Y[1,5]^3
        Y[14,18] = - A*Y[14,5]/Y[1,5]^2 + 2.0*A*Y[3,5]*Y[5,5]/Y[1,5]^3
        Y[16,18] = Y[16,4]- A*Y[16,5]/Y[1,5]^2 + 2.0*A*Y[5,5]^2 / Y[1,5]^3

        # Function Y19
        
        A = -3.496e+3
        Y[1,19]  = A*C[1,19]
        Y[2,19]  = A*C[2,19]
        Y[3,19]  = A*C[3,19]
        Y[5,19]  = A*C[5,19]
        Y[7,19]  = A*C[7,19]
        Y[8,19]  = A*C[8,19]
        Y[9,19]  = A*C[9,19]
        Y[13,19] = A*C[13,19]
        Y[14,19] = A*C[14,19]
        Y[16,19] = A*C[16,19]

        # Function Y20
        
        A = 6.2212e+4
        Y[1,20]  = A/C[1,17] - Y[1,1]
        Y[2,20]  = - A*C[2,17]/C[1,17]^2
        Y[3,20]  = - A*C[3,17]/C[1,17]^2 - Y[3,1]
        Y[4,20]  = - A*C[4,17]/C[1,17]^2 - Y[4,1]
        Y[5,20]  = - A*C[5,17]/C[1,17]^2
        Y[6,20]  = - A*C[6,17]/C[1,17]^2
        Y[7,20]  = - A*C[7,17]/C[1,17]^2 + 2.0*A*C[2,17]^2/C[1,17]^3
        Y[8,20]  = - A*C[8,17]/C[1,17]^2 + 2.0*A*C[2,17]*C[3,17]/C[1,17]^3
        Y[9,20]  = - A*C[9,17]/C[1,17]^2 + 2.0*A*C[3,17]^2/C[1,17]^3
        Y[10,20] = - A*C[10,17]/C[1,17]^2 + 2.0*A*C[2,17]*C[4,17]/C[1,17]^3
        Y[11,20] = - A*C[11,17]/C[1,17]^2 + 2.0*A*C[3,17]*C[4,17]/C[1,17]^3
        Y[12,20] = - A*C[12,17]/C[1,17]^2 + 2.0*A*C[4,17]^2/C[1,17]^3
        Y[13,20] = - A*C[13,17]/C[1,17]^2 + 2.0*A*C[2,17]*C[5,17]/C[1,17]^3
        Y[14,20] = - A*C[14,17]/C[1,17]^2 + 2.0*A*C[3,17]*C[5,17]/C[1,17]^3
        Y[15,20] = - A*C[15,17]/C[1,17]^2 + 2.0*A*C[4,17]*C[5,17]/C[1,17]^3
        Y[16,20] = - A*C[16,17]/C[1,17]^2 + 2.0*A*C[5,17]^2/C[1,17]^3
        Y[17,20] = - A*C[17,17]/C[1,17]^2 + 2.0*A*C[2,17]*C[6,17]/C[1,17]^3
        Y[18,20] = - A*C[18,17]/C[1,17]^2 + 2.0*A*C[3,17]*C[6,17]/C[1,17]^3
        Y[19,20] = - A*C[19,17]/C[1,17]^2 + 2.0*A*C[4,17]*C[6,17]/C[1,17]^3
        Y[20,20] = - A*C[20,17]/C[1,17]^2 + 2.0*A*C[5,17]*C[6,17]/C[1,17]^3
        Y[21,20] = - A*C[21,17]/C[1,17]^2 + 2.0*A*C[6,17]^2/C[1,17]^3

        return C, Y

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

